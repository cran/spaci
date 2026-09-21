#' Rescale a numeric object to the unit interval
#'
#' Linearly rescales a numeric vector or matrix to `[0, 1]`. If the input has no
#' finite range (constant or all non-finite), a zero object of the same shape is
#' returned so that a degenerate distance component contributes nothing to a
#' composite distance.
#'
#' @param x A numeric vector or matrix.
#' @return A numeric object of the same shape as `x`, rescaled to `[0, 1]`.
#' @keywords internal
#' @noRd
range01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || isTRUE(all.equal(rng[1], rng[2]))) {
    if (is.matrix(x)) {
      return(matrix(0, nrow = nrow(x), ncol = ncol(x)))
    }
    return(rep(0, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

#' Clip propensity scores away from 0 and 1
#'
#' Non-finite values are replaced by 0.5 and the result is clipped to
#' `[eps, 1 - eps]` to keep inverse-probability weights finite.
#'
#' @param p Numeric vector of probabilities.
#' @param eps Small positive truncation constant.
#' @return Numeric vector clipped to `[eps, 1 - eps]`.
#' @keywords internal
#' @noRd
clip_ps <- function(p, eps = 1e-6) {
  p <- as.numeric(p)
  p[!is.finite(p)] <- 0.5
  pmin(pmax(p, eps), 1 - eps)
}

#' Standardise a numeric vector, robust to degenerate input
#'
#' Centres and scales `x` to unit standard deviation. If the standard deviation
#' is zero or non-finite, a zero vector is returned.
#'
#' @param x Numeric vector.
#' @return Standardised numeric vector.
#' @keywords internal
#' @noRd
safe_scale <- function(x) {
  x <- as.numeric(x)
  sx <- stats::sd(x, na.rm = TRUE)
  mx <- mean(x, na.rm = TRUE)
  if (!is.finite(sx) || sx == 0) {
    return(rep(0, length(x)))
  }
  out <- (x - mx) / sx
  out[!is.finite(out)] <- 0
  out
}

#' Validate and coerce the common estimator inputs
#'
#' Shared input handling for the estimators. Coerces the covariate argument to a
#' numeric matrix with column names and checks that the treatment is binary and
#' the dimensions are consistent.
#'
#' @param Y Numeric outcome vector.
#' @param Z Binary treatment vector (0/1).
#' @param X Covariate matrix or data frame (may be `NULL`).
#' @param coords Two-column matrix or data frame of spatial coordinates
#'   (may be `NULL` when coordinates are not required).
#' @param require_X,require_coords Logical; error if the corresponding argument
#'   is missing when it is required by the caller.
#' @return A list with cleaned `Y`, `Z`, `X` (matrix or `NULL`) and `coords`.
#' @keywords internal
#' @noRd
validate_inputs <- function(Y, Z, X = NULL, coords = NULL,
                            require_X = FALSE, require_coords = FALSE) {
  if (!is.numeric(Y) && !is.logical(Y)) {
    stop("`Y` must be a numeric vector.", call. = FALSE)
  }
  Y <- as.numeric(Y)
  Z <- as.numeric(Z)
  n <- length(Y)

  if (n < 4L) {
    stop("At least four observations are required.", call. = FALSE)
  }
  if (length(Z) != n) {
    stop("`Y` and `Z` must have the same length.", call. = FALSE)
  }
  if (anyNA(Y) || anyNA(Z)) {
    stop("`Y` and `Z` must not contain missing values.", call. = FALSE)
  }
  if (!all(is.finite(Y))) {
    stop("`Y` must contain only finite values.", call. = FALSE)
  }
  if (!all(Z %in% c(0, 1))) {
    stop("`Z` must be a binary (0/1) treatment vector.", call. = FALSE)
  }
  if (length(unique(Z)) < 2) {
    stop("`Z` must contain both treated and control units.", call. = FALSE)
  }
  if (min(sum(Z == 1), sum(Z == 0)) < 2L) {
    stop("At least two treated and two control units are required.",
         call. = FALSE)
  }

  if (require_X && is.null(X)) {
    stop("`X` (covariates) must be supplied.", call. = FALSE)
  }
  if (!is.null(X)) {
    if (is.data.frame(X)) {
      X <- data.matrix(X)
    } else if (is.vector(X)) {
      X <- matrix(X, ncol = 1L)
    }
    X <- as.matrix(X)
    if (nrow(X) != n) {
      stop("`X` must have one row per observation.", call. = FALSE)
    }
    if (anyNA(X) || !all(is.finite(X))) {
      stop("`X` must not contain missing or non-finite values.", call. = FALSE)
    }
    if (is.null(colnames(X))) {
      colnames(X) <- paste0("X", seq_len(ncol(X)))
    }
    storage.mode(X) <- "double"
  }

  if (require_coords && is.null(coords)) {
    stop("`coords` (spatial coordinates) must be supplied.", call. = FALSE)
  }
  if (!is.null(coords)) {
    coords <- as.matrix(coords)
    if (nrow(coords) != n) {
      stop("`coords` must have one row per observation.", call. = FALSE)
    }
    if (ncol(coords) != 2) {
      stop("`coords` must have exactly two columns (x, y).", call. = FALSE)
    }
    if (anyNA(coords) || !all(is.finite(coords))) {
      stop("`coords` must not contain missing or non-finite values.",
           call. = FALSE)
    }
    storage.mode(coords) <- "double"
  }

  list(Y = Y, Z = Z, X = X, coords = coords)
}

#' Argument checks shared by the estimators
#'
#' @param caliper Matching caliper (positive, may be `Inf`).
#' @param tau Kernel bandwidth (positive finite), or `NULL` to skip.
#' @param level Confidence level in (0, 1).
#' @param grid_step Grid step in (0, 1], or `NULL` to skip.
#' @return `NULL`, invisibly; called for its side effect of erroring.
#' @keywords internal
#' @noRd
check_scalars <- function(caliper = NULL, tau = NULL, level = NULL,
                          grid_step = NULL) {
  if (!is.null(caliper) &&
      (!is.numeric(caliper) || length(caliper) != 1L || is.na(caliper) ||
       caliper <= 0)) {
    stop("`caliper` must be a single positive number (or Inf).", call. = FALSE)
  }
  if (!is.null(tau) &&
      (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0)) {
    stop("`tau` must be a single positive number.", call. = FALSE)
  }
  if (!is.null(level) &&
      (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
       level <= 0 || level >= 1)) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.null(grid_step) &&
      (!is.numeric(grid_step) || length(grid_step) != 1L ||
       !is.finite(grid_step) || grid_step <= 0 || grid_step > 1)) {
    stop("grid step must be a single number in (0, 1].", call. = FALSE)
  }
  invisible(NULL)
}
