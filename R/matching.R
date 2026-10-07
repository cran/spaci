#' Greedy 1:1 nearest-neighbour matching and ATT
#'
#' Matches each treated unit to its nearest available control under a supplied
#' distance matrix (without replacement, within an optional caliper) and returns
#' the average treatment effect on the treated as the mean of within-pair
#' outcome differences.
#'
#' Treated units are processed in random order, so the result depends on the RNG
#' state; set a seed for reproducibility.
#'
#' @param Dmat Square distance matrix between all units.
#' @param Y Numeric outcome vector.
#' @param Z Binary treatment vector (0/1).
#' @param caliper Maximum acceptable matching distance. Pairs whose distance
#'   exceeds the caliper are discarded. Defaults to `Inf` (no caliper).
#'
#' @return A list with `att`, `se`, the number of matched (`n_match`) and
#'   dropped (`n_drop`) treated units, and the matched `pairs` (a two-column
#'   matrix of treated/control row indices).
#' @param method Matching algorithm: `"greedy"` (default) processes treated
#'   units in random order, matching each to its nearest available control
#'   (depends on the RNG); `"optimal"` solves the 1:1 assignment minimising total
#'   distance via [clue::solve_LSAP()] (deterministic, requires the `clue`
#'   package).
#' @keywords internal
#' @noRd
match_att <- function(Dmat, Y, Z, caliper = Inf, method = c("greedy", "optimal")) {
  method <- match.arg(method)
  treated <- which(Z == 1)
  controls <- which(Z == 0)

  empty_fit <- function(n_match = 0, n_drop = length(treated), pairs = NULL) {
    list(att = NA_real_, se = NA_real_, n_match = n_match,
         n_drop = n_drop, pairs = pairs)
  }

  if (length(treated) == 0 || length(controls) == 0) {
    return(empty_fit())
  }

  pairs <- if (method == "optimal") {
    match_pairs_optimal(Dmat, treated, controls, caliper)
  } else {
    match_pairs_greedy(Dmat, treated, controls, caliper)
  }

  if (is.null(pairs) || nrow(pairs) < 2) {
    return(empty_fit(n_match = if (is.null(pairs)) 0 else nrow(pairs),
                     n_drop = length(treated) -
                       (if (is.null(pairs)) 0 else nrow(pairs)),
                     pairs = pairs))
  }

  diffs <- Y[pairs[, 1]] - Y[pairs[, 2]]
  list(att = mean(diffs), se = stats::sd(diffs) / sqrt(length(diffs)),
       n_match = nrow(pairs), n_drop = length(treated) - nrow(pairs),
       pairs = pairs)
}

#' Greedy 1:1 nearest-neighbour matching (random treated order)
#' @keywords internal
#' @noRd
match_pairs_greedy <- function(Dmat, treated, controls, caliper) {
  available_controls <- controls
  pairs <- matrix(integer(0), ncol = 2)
  treated_order <- sample(treated, length(treated), replace = FALSE)

  for (i in treated_order) {
    if (length(available_controls) == 0) break
    dvec <- Dmat[i, available_controls]
    if (all(!is.finite(dvec))) next
    jpos <- which.min(dvec)
    dmin <- dvec[jpos]
    if (!is.finite(dmin) || dmin > caliper) next
    j <- available_controls[jpos]
    pairs <- rbind(pairs, c(i, j))
    available_controls <- setdiff(available_controls, j)
  }
  pairs
}

#' Optimal 1:1 assignment matching (deterministic, via clue::solve_LSAP)
#' @keywords internal
#' @noRd
match_pairs_optimal <- function(Dmat, treated, controls, caliper) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("`match_method = \"optimal\"` requires the 'clue' package.", call. = FALSE)
  }
  ## assign the smaller group to the larger so solve_LSAP has nrow <= ncol
  swap <- length(treated) > length(controls)
  rows <- if (swap) controls else treated
  cols <- if (swap) treated else controls

  cost <- Dmat[rows, cols, drop = FALSE]
  big <- max(cost[is.finite(cost)], na.rm = TRUE) * length(cost) + 1
  cost[!is.finite(cost)] <- big

  assign <- tryCatch(as.integer(clue::solve_LSAP(cost)),
                     error = function(e) NULL)
  if (is.null(assign)) return(matrix(integer(0), ncol = 2))

  ri <- seq_along(rows)
  d_sel <- cost[cbind(ri, assign)]
  keep <- is.finite(d_sel) & d_sel <= caliper & d_sel < big
  t_idx <- if (swap) cols[assign[keep]] else rows[ri[keep]]
  c_idx <- if (swap) rows[ri[keep]] else cols[assign[keep]]
  cbind(t_idx, c_idx)
}

#' Absolute standardised mean difference in covariates over matched pairs
#'
#' @param pairs Two-column matrix of matched treated/control indices.
#' @param X Numeric covariate matrix.
#' @return The summed absolute standardised mean difference, or `Inf` if fewer
#'   than two pairs are available.
#' @keywords internal
#' @noRd
covariate_balance_score <- function(pairs, X) {
  if (is.null(pairs) || nrow(pairs) < 2) return(Inf)

  treated_idx <- pairs[, 1]
  control_idx <- pairs[, 2]

  X_T <- X[treated_idx, , drop = FALSE]
  X_C <- X[control_idx, , drop = FALSE]

  sX <- apply(X_T, 2, stats::sd, na.rm = TRUE)
  sX[!is.finite(sX) | sX == 0] <- 1

  sum(abs((colMeans(X_T) - colMeans(X_C)) / sX))
}

#' Composite balance score used to select the iDAPS weights
#'
#' Combines covariate balance, mean spatial proximity of matched pairs, and
#' balance in neighbourhood exposure, as in the objective of Equation (2.5).
#'
#' @param pairs Two-column matrix of matched treated/control indices.
#' @param X Numeric covariate matrix.
#' @param E Numeric vector of neighbourhood exposures.
#' @param d_space_raw Raw (unnormalised) Euclidean distance matrix.
#' @return The combined balance score (lower is better), or `Inf` if fewer than
#'   two pairs are available.
#' @keywords internal
#' @noRd
balance_score_idaps <- function(pairs, X, E, d_space_raw) {
  if (is.null(pairs) || nrow(pairs) < 2) return(Inf)

  treated_idx <- pairs[, 1]
  control_idx <- pairs[, 2]

  BX <- covariate_balance_score(pairs, X)

  dbar <- mean(d_space_raw[cbind(treated_idx, control_idx)], na.rm = TRUE)
  if (!is.finite(dbar)) dbar <- Inf

  sG <- stats::sd(E[treated_idx], na.rm = TRUE)
  if (!is.finite(sG) || sG == 0) sG <- 1
  BG <- abs((mean(E[treated_idx]) - mean(E[control_idx])) / sG)

  BX + dbar + BG
}
