#' Naive propensity-score matching ATT
#'
#' Estimates the ATT by greedy 1:1 matching on the (normalised) difference in
#' estimated propensity scores. This is the baseline that ignores both spatial
#' confounding and spatial interference.
#'
#' @param Y Numeric outcome vector.
#' @param Z Binary treatment vector (0/1).
#' @param X Covariate matrix or data frame.
#' @param caliper Maximum acceptable matching distance (default `0.25`).
#' @param match_method Matching algorithm: `"greedy"` (default) matches treated
#'   units in random order to their nearest available control (depends on
#'   `seed`); `"optimal"` solves the 1:1 assignment minimising total distance
#'   (deterministic, needs the `clue` package).
#' @param seed Optional integer seed for the (randomised) matching order.
#' @param level Confidence level for the reported interval.
#'
#' @return An `idaps_fit` object.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 1)
#' naive_ps(sim$Y, sim$Z, sim$X)
#' @export
naive_ps <- function(Y, Z, X, caliper = 0.25, match_method = c("greedy", "optimal"),
                     seed = NULL, level = 0.95) {
  match_method <- match.arg(match_method)
  check_scalars(caliper = caliper, level = level)
  if (!is.null(seed)) {
    old <- .save_seed(); on.exit(.restore_seed(old), add = TRUE); set.seed(seed)
  }
  d <- validate_inputs(Y, Z, X, require_X = TRUE)
  ps <- estimate_ps(d$Z, d$X)

  D_ps <- range01(abs(outer(ps, ps, "-")))
  diag(D_ps) <- 0

  fit <- match_att(D_ps, d$Y, d$Z, caliper, method = match_method)
  new_idaps_fit("Naive PS", fit$att, fit$se, fit$n_match, fit$n_drop,
                extras = list(ps = ps, pairs = fit$pairs), level = level)
}

#' Distance-adjusted propensity score (DAPS) matching ATT
#'
#' Estimates the ATT by matching on a convex combination of propensity-score
#' distance and spatial (Euclidean) distance,
#' \eqn{\alpha D^{PS} + (1 - \alpha) D^{Spatial}}, adjusting for spatial
#' confounding but not interference. The weight \eqn{\alpha} is chosen over a
#' grid to minimise covariate imbalance in the matched sample
#' (Papadogeorgou et al., 2019).
#'
#' @inheritParams naive_ps
#' @param coords Two-column matrix or data frame of spatial coordinates.
#' @param alpha_grid_step Grid step for the \eqn{\alpha} search over `[0, 1]`.
#'
#' @return An `idaps_fit` object; the selected \eqn{\alpha} is
#'   reported in `weights`.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 1)
#' daps(sim$Y, sim$Z, sim$X, sim$coords)
#' @export
daps <- function(Y, Z, X, coords, caliper = 0.25, alpha_grid_step = 0.1,
                 match_method = c("greedy", "optimal"), seed = NULL,
                 level = 0.95) {
  match_method <- match.arg(match_method)
  check_scalars(caliper = caliper, level = level, grid_step = alpha_grid_step)
  if (!is.null(seed)) {
    old <- .save_seed(); on.exit(.restore_seed(old), add = TRUE); set.seed(seed)
  }
  d <- validate_inputs(Y, Z, X, coords, require_X = TRUE, require_coords = TRUE)
  ps <- estimate_ps(d$Z, d$X)
  comp <- distance_components(ps, d$coords, rep(0, length(d$Y)))

  grid_vals <- seq(0, 1, by = alpha_grid_step)
  best_B <- Inf; best_alpha <- NA_real_; best_fit <- NULL

  for (alpha in grid_vals) {
    D <- alpha * comp$D_ps + (1 - alpha) * comp$d_space
    diag(D) <- 0
    fit <- match_att(D, d$Y, d$Z, caliper, method = match_method)
    if (is.na(fit$att) || is.null(fit$pairs) || nrow(fit$pairs) < 2) next
    B <- covariate_balance_score(fit$pairs, d$X)
    if (is.finite(B) && B < best_B) {
      best_B <- B; best_alpha <- alpha; best_fit <- fit
    }
  }

  if (is.null(best_fit)) {
    return(new_idaps_fit("DAPS", NA_real_, NA_real_, 0, sum(d$Z == 1),
                         weights = c(alpha = NA_real_), level = level))
  }

  new_idaps_fit("DAPS", best_fit$att, best_fit$se, best_fit$n_match,
                best_fit$n_drop, weights = c(alpha = best_alpha),
                extras = list(ps = ps, pairs = best_fit$pairs,
                              balance = best_B), level = level)
}

#' Distance-adjusted propensity score with interference (iDAPS)
#'
#' Estimates the ATT by matching on the composite distance of Equation (2.3),
#' \eqn{\pi_1 D^{PS} + \pi_2 D^{Spatial} + \pi_3 D^{Interference}}, where the
#' interference component is the absolute difference in neighbourhood exposure.
#' The weights \eqn{(\pi_1, \pi_2, \pi_3)} sum to one and are chosen over a grid
#' to minimise a composite balance score (covariate balance, spatial proximity
#' and exposure balance), rather than tuned by hand. Setting \eqn{\pi_3 = 0}
#' recovers DAPS and \eqn{\pi_2 = \pi_3 = 0} recovers naive PS.
#'
#' @inheritParams daps
#' @param tau Positive spatial decay parameter of the exposure kernel (default
#'   `0.1`); see [neighbourhood_exposure()].
#' @param normalize Logical; row-normalise the exposure kernel (default `TRUE`).
#' @param pi_grid_step Grid step for the weight search over the simplex.
#'
#' @return An `idaps_fit` object; the selected weights are
#'   reported in `weights`.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 1)
#' idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1)
#' @export
idaps <- function(Y, Z, X, coords, tau = 0.1, caliper = 0.25,
                  pi_grid_step = 0.1, normalize = TRUE,
                  match_method = c("greedy", "optimal"), seed = NULL,
                  level = 0.95) {
  match_method <- match.arg(match_method)
  check_scalars(caliper = caliper, tau = tau, level = level,
                grid_step = pi_grid_step)
  if (!isTRUE(normalize) && !isFALSE(normalize)) {
    stop("`normalize` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(seed)) {
    old <- .save_seed(); on.exit(.restore_seed(old), add = TRUE); set.seed(seed)
  }
  d <- validate_inputs(Y, Z, X, coords, require_X = TRUE, require_coords = TRUE)
  ps <- estimate_ps(d$Z, d$X)
  E <- neighbourhood_exposure(d$coords, d$Z, tau = tau, normalize = normalize)$E
  comp <- distance_components(ps, d$coords, E)

  grid_vals <- seq(0, 1, by = pi_grid_step)
  best_B <- Inf; best_pi <- c(NA_real_, NA_real_, NA_real_); best_fit <- NULL

  for (pi1 in grid_vals) {
    for (pi2 in grid_vals) {
      pi3 <- 1 - pi1 - pi2
      if (pi3 < -1e-12) next
      if (pi3 < 0) pi3 <- 0

      D <- pi1 * comp$D_ps + pi2 * comp$d_space + pi3 * comp$d_exp
      diag(D) <- 0
      fit <- match_att(D, d$Y, d$Z, caliper, method = match_method)
      if (is.na(fit$att) || is.null(fit$pairs) || nrow(fit$pairs) < 2) next
      B <- balance_score_idaps(fit$pairs, d$X, E, comp$d_space_raw)
      if (is.finite(B) && B < best_B) {
        best_B <- B; best_pi <- c(pi1, pi2, pi3); best_fit <- fit
      }
    }
  }

  if (is.null(best_fit)) {
    return(new_idaps_fit("iDAPS", NA_real_, NA_real_, 0, sum(d$Z == 1),
                         weights = c(pi1 = NA_real_, pi2 = NA_real_,
                                     pi3 = NA_real_), level = level))
  }

  new_idaps_fit("iDAPS", best_fit$att, best_fit$se, best_fit$n_match,
                best_fit$n_drop,
                weights = c(pi1 = best_pi[1], pi2 = best_pi[2],
                            pi3 = best_pi[3]),
                extras = list(ps = ps, E = E, pairs = best_fit$pairs,
                              balance = best_B), level = level)
}

## --- RNG state helpers -------------------------------------------------------
.save_seed <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}
.restore_seed <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", old, envir = globalenv())
  }
}
