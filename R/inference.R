#' Data-driven bandwidth for the spatial HAC variance
#'
#' Estimates the spatial dependence range of the influence values by binning the
#' pairwise products by distance and finding where the autocorrelation decays
#' below a threshold. Falls back to one tenth of the maximum pairwise distance.
#'
#' @param infl Numeric influence values.
#' @param coords Two-column coordinate matrix.
#' @param threshold Autocorrelation cutoff (default `0.1`).
#' @return A single positive bandwidth.
#' @keywords internal
#' @noRd
auto_bandwidth <- function(infl, coords, threshold = 0.1) {
  z <- infl - mean(infl)
  v <- mean(z^2)
  D <- as.matrix(stats::dist(coords))
  maxd <- max(D)
  fallback <- maxd / 10
  if (!is.finite(v) || v <= 0) return(fallback)

  ut <- upper.tri(D)
  dvec <- D[ut]
  pvec <- outer(z, z)[ut] / v
  brks <- seq(0, maxd, length.out = 13)
  mids <- (brks[-1] + brks[-length(brks)]) / 2
  ac <- vapply(seq_along(mids), function(k) {
    sel <- dvec > brks[k] & dvec <= brks[k + 1]
    if (sum(sel) < 5) NA_real_ else mean(pvec[sel])
  }, numeric(1))

  below <- which(is.finite(ac) & ac < threshold)
  if (length(below) == 0) return(fallback)
  mids[below[1]]
}

#' Spatial HAC (Conley) variance for a doubly robust ATT
#'
#' Computes a spatial heteroskedasticity- and autocorrelation-consistent
#' (Conley) variance for a [recoverU()] or [recoverUplus()] estimate,
#' \deqn{\widehat{\mathrm{Var}}(\hat\tau) = \frac{1}{n^2} \sum_i \sum_j
#'   K\!\left(\frac{\|s_i - s_j\|}{b}\right) \hat\psi_i \hat\psi_j,}
#' where \eqn{\hat\psi_i} are the estimator's influence values and \eqn{K} is a
#' kernel with bandwidth \eqn{b}. Unlike the default i.i.d. standard error, this
#' accounts for spatial correlation of the influence values, which otherwise
#' makes intervals too narrow.
#'
#' @param object A `recoverU`/`recoverU+` `idaps_fit` (carries the influence
#'   values and coordinates). Matching estimators have no influence-value
#'   representation; use [boot_spatial()] for those.
#' @param bandwidth Kernel bandwidth: `"auto"` (default, estimated from the
#'   influence-value autocorrelation range) or a positive number.
#' @param kernel `"bartlett"` (default, guarantees a non-negative variance) or
#'   `"uniform"`.
#' @param level Confidence level for the returned interval.
#'
#' @return An object of class `spaci_vcov`: a list with the HAC `variance`, `se`,
#'   confidence interval `ci`, the `att`, and the `bandwidth`/`kernel` used.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 1)
#' fit <- recoverUplus(sim$Y, sim$Z, sim$X, sim$coords)
#' vcov_hac(fit)
#' @export
vcov_hac <- function(object, bandwidth = "auto",
                     kernel = c("bartlett", "uniform"), level = 0.95) {
  kernel <- match.arg(kernel)
  check_scalars(level = level)
  if (!inherits(object, "idaps_fit")) {
    stop("`object` must be an idaps_fit.", call. = FALSE)
  }
  infl <- object$extras$infl
  coords <- object$extras$coords
  if (is.null(infl) || is.null(coords)) {
    stop("No influence values in `object`; vcov_hac() applies to recoverU / ",
         "recoverU+ fits. Use boot_spatial() for matching estimators.",
         call. = FALSE)
  }

  if (identical(bandwidth, "auto")) {
    bandwidth <- auto_bandwidth(infl, coords)
  }
  if (!is.numeric(bandwidth) || length(bandwidth) != 1L || bandwidth <= 0) {
    stop("`bandwidth` must be \"auto\" or a single positive number.",
         call. = FALSE)
  }

  n <- length(infl)
  D <- as.matrix(stats::dist(coords))
  u <- D / bandwidth
  K <- switch(kernel,
              bartlett = pmax(0, 1 - u),
              uniform  = (u <= 1) * 1)
  variance <- sum(K * outer(infl, infl)) / n^2
  variance <- max(variance, 0)
  se <- sqrt(variance)
  zq <- stats::qnorm(1 - (1 - level) / 2)

  structure(
    list(method = object$method, att = object$att, variance = variance,
         se = se, ci = c(lower = object$att - zq * se,
                         upper = object$att + zq * se),
         level = level, bandwidth = bandwidth, kernel = kernel,
         se_iid = object$se),
    class = "spaci_vcov")
}

#' Print method for a spatial HAC variance
#' @param x A `spaci_vcov` object.
#' @param digits Significant digits.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.spaci_vcov <- function(x, digits = 4, ...) {
  cat("Spatial HAC (Conley) inference for", x$method, "\n\n")
  cat("  ATT     =", format(x$att, digits = digits), "\n")
  cat("  SE(HAC) =", format(x$se, digits = digits),
      sprintf("(i.i.d. SE = %s)\n", format(x$se_iid, digits = digits)))
  cat(sprintf("  %g%% CI  = [%s, %s]\n", 100 * x$level,
              format(x$ci[["lower"]], digits = digits),
              format(x$ci[["upper"]], digits = digits)))
  cat(sprintf("  kernel = %s, bandwidth = %s\n", x$kernel,
              format(x$bandwidth, digits = digits)))
  invisible(x)
}

## -- internal estimator dispatch used by the resampling procedures ------------
run_estimator <- function(method, Y, Z, X, coords, tau, caliper,
                          matern_method, matern_nu, seed = NULL) {
  switch(method,
    naive_ps     = naive_ps(Y, Z, X, caliper = caliper, seed = seed),
    daps         = daps(Y, Z, X, coords, caliper = caliper, seed = seed),
    idaps        = idaps(Y, Z, X, coords, tau = tau, caliper = caliper, seed = seed),
    recoverU     = recoverU(Y, Z, X, coords, tau = tau,
                            matern_method = matern_method, matern_nu = matern_nu),
    recoverUplus = recoverUplus(Y, Z, X, coords, tau = tau,
                                matern_method = matern_method, matern_nu = matern_nu),
    stop("Unknown `method`: ", method, call. = FALSE))
}

#' Spatial block bootstrap confidence interval
#'
#' Resamples square spatial blocks of units with replacement to build a
#' bootstrap distribution that reflects spatial dependence. For the doubly robust
#' estimators (`recoverU`/`recoverU+`) it block-resamples the estimator's
#' influence values (a fast, stable linearised bootstrap that avoids re-kriging
#' resampled coordinates); for the matching estimators it re-runs the whole
#' pipeline on each resample (recomputing neighbourhood exposures, with a tiny
#' coordinate jitter to avoid exact-duplicate ties).
#'
#' @param Y,Z,X,coords Data, as passed to the estimators.
#' @param method Estimator to bootstrap: one of `"recoverUplus"`, `"recoverU"`,
#'   `"idaps"`, `"daps"`, `"naive_ps"`.
#' @param block_size Side length of the square blocks: `"auto"` (default,
#'   `domain width x n^(-1/4)`, giving on the order of sqrt(n) blocks) or a
#'   positive number. The block side should exceed the spatial dependence range.
#' @param B Number of bootstrap resamples (default `200`).
#' @param tau,caliper,matern_method,matern_nu Estimator arguments (as relevant).
#' @param level Confidence level.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class `spaci_boot` with the bootstrap `se`, percentile
#'   `ci`, point estimate, and the retained resample estimates.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 1)
#' boot_spatial(sim$Y, sim$Z, sim$X, sim$coords, method = "recoverUplus", B = 50)
#' @export
boot_spatial <- function(Y, Z, X, coords,
                         method = c("recoverUplus", "recoverU", "idaps",
                                    "daps", "naive_ps"),
                         block_size = "auto", B = 200, tau = 0.1,
                         caliper = 0.25, matern_method = c("mle", "geoR"),
                         matern_nu = 0.5, level = 0.95, seed = NULL) {
  method <- match.arg(method)
  matern_method <- match.arg(matern_method)
  check_scalars(caliper = caliper, tau = tau, level = level)
  if (!is.null(seed)) {
    old <- .save_seed(); on.exit(.restore_seed(old), add = TRUE); set.seed(seed)
  }
  d <- validate_inputs(Y, Z, X, coords, require_X = TRUE, require_coords = TRUE)
  Y <- d$Y; Z <- d$Z; X <- d$X; coords <- d$coords
  n <- length(Y)

  dom <- max(apply(coords, 2, function(x) diff(range(x))))
  if (identical(block_size, "auto")) {
    ## standard growing-block rule: side proportional to n^(-1/4) gives on the
    ## order of sqrt(n) blocks, with side -> Inf and side/domain -> 0, and large
    ## enough to exceed the spatial dependence range in typical designs.
    block_size <- dom * length(Y)^(-0.25)
  }
  if (!is.numeric(block_size) || block_size <= 0) {
    stop("`block_size` must be \"auto\" or a positive number.", call. = FALSE)
  }
  jit <- block_size / 50

  ## block id per unit; list of unit indices per block
  bx <- floor((coords[, 1] - min(coords[, 1])) / block_size)
  by <- floor((coords[, 2] - min(coords[, 2])) / block_size)
  bid <- paste(bx, by)
  blocks <- split(seq_len(n), bid)
  nb <- length(blocks)

  point_fit <- run_estimator(method, Y, Z, X, coords, tau, caliper,
                             matern_method, matern_nu, seed = 1)
  point <- point_fit$att
  infl <- point_fit$extras$infl          # non-NULL for recoverU / recoverU+

  if (!is.null(infl)) {
    ## linearised block bootstrap: resample blocks of (centred) influence values
    infl <- infl - mean(infl)
    ests <- vapply(seq_len(B), function(b) {
      idx <- unlist(blocks[sample(nb, nb, replace = TRUE)], use.names = FALSE)
      point + mean(infl[idx])
    }, numeric(1))
  } else {
    ## re-run block bootstrap for matching estimators (no kriging step)
    ests <- rep(NA_real_, B)
    for (b in seq_len(B)) {
      idx <- unlist(blocks[sample(nb, nb, replace = TRUE)], use.names = FALSE)
      cc <- coords[idx, , drop = FALSE] +
        matrix(stats::rnorm(length(idx) * 2, 0, jit), ncol = 2)
      fit <- tryCatch(
        run_estimator(method, Y[idx], Z[idx], X[idx, , drop = FALSE], cc,
                      tau, caliper, matern_method, matern_nu, seed = 1),
        error = function(e) NULL)
      if (!is.null(fit) && is.finite(fit$att)) ests[b] <- fit$att
    }
  }
  ests <- ests[is.finite(ests)]
  a <- (1 - level) / 2
  ci <- stats::quantile(ests, c(a, 1 - a), names = FALSE)

  structure(
    list(method = method, att = point, se = stats::sd(ests),
         ci = c(lower = ci[1], upper = ci[2]), level = level,
         B = length(ests), block_size = block_size, estimates = ests),
    class = "spaci_boot")
}

#' Print method for a spatial block bootstrap
#' @param x A `spaci_boot` object.
#' @param digits Significant digits.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.spaci_boot <- function(x, digits = 4, ...) {
  cat("Spatial block bootstrap for", x$method,
      sprintf("(%d resamples)\n\n", x$B))
  cat("  ATT       =", format(x$att, digits = digits), "\n")
  cat("  SE(boot)  =", format(x$se, digits = digits), "\n")
  cat(sprintf("  %g%% CI    = [%s, %s]  (percentile)\n", 100 * x$level,
              format(x$ci[["lower"]], digits = digits),
              format(x$ci[["upper"]], digits = digits)))
  cat("  block size =", format(x$block_size, digits = digits), "\n")
  invisible(x)
}

#' Conditional randomization test for no direct effect
#'
#' Design-based test of the sharp null that the treatment has no direct effect on
#' any unit. Treatment is redrawn from the estimated propensity score (so the
#' reference distribution respects the observed treatment mechanism),
#' neighbourhood exposures are recomputed, and the estimator is recomputed for
#' each draw. The two-sided p-value is
#' \eqn{(1 + \#\{|\hat\tau^{(r)}| \ge |\hat\tau|\}) / (1 + R)}.
#'
#' @inheritParams boot_spatial
#' @param R Number of randomization draws (default `200`).
#'
#' @return An object of class `spaci_randtest` with the observed statistic, the
#'   p-value and the null draws.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 1)
#' rand_test(sim$Y, sim$Z, sim$X, sim$coords, method = "naive_ps", R = 50)
#' \donttest{
#' # refits the full recoverU+ pipeline on every draw, so it is slower
#' rand_test(sim$Y, sim$Z, sim$X, sim$coords, method = "recoverUplus", R = 20)
#' }
#' @export
rand_test <- function(Y, Z, X, coords,
                      method = c("recoverUplus", "recoverU", "idaps",
                                 "daps", "naive_ps"),
                      R = 200, tau = 0.1, caliper = 0.25,
                      matern_method = c("mle", "geoR"), matern_nu = 0.5,
                      seed = NULL) {
  method <- match.arg(method)
  matern_method <- match.arg(matern_method)
  check_scalars(caliper = caliper, tau = tau)
  if (!is.null(seed)) {
    old <- .save_seed(); on.exit(.restore_seed(old), add = TRUE); set.seed(seed)
  }
  d <- validate_inputs(Y, Z, X, coords, require_X = TRUE, require_coords = TRUE)
  Y <- d$Y; Z <- d$Z; X <- d$X; coords <- d$coords

  ps <- estimate_ps(Z, X)                       # propensity given covariates
  obs <- abs(run_estimator(method, Y, Z, X, coords, tau, caliper,
                           matern_method, matern_nu, seed = 1)$att)

  null <- rep(NA_real_, R)
  for (r in seq_len(R)) {
    Zr <- stats::rbinom(length(Z), 1, ps)
    if (length(unique(Zr)) < 2) next
    fit <- tryCatch(
      run_estimator(method, Y, Zr, X, coords, tau, caliper,
                    matern_method, matern_nu, seed = 1),
      error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$att)) null[r] <- abs(fit$att)
  }
  null <- null[is.finite(null)]
  pval <- (1 + sum(null >= obs)) / (1 + length(null))

  structure(list(method = method, statistic = obs, p_value = pval,
                 R = length(null), null = null),
            class = "spaci_randtest")
}

#' Print method for a randomization test
#' @param x A `spaci_randtest` object.
#' @param digits Significant digits.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.spaci_randtest <- function(x, digits = 4, ...) {
  cat("Conditional randomization test (H0: no direct effect) for", x$method, "\n\n")
  cat("  |ATT observed| =", format(x$statistic, digits = digits), "\n")
  cat("  p-value        =", format(x$p_value, digits = digits),
      sprintf("(%d draws)\n", x$R))
  invisible(x)
}
