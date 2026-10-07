#' Matern covariance matrix
#'
#' Evaluates a Matern covariance matrix for a matrix of distances, using the
#' parametrisation \eqn{C(h) = \sigma^2 \frac{(z)^\nu K_\nu(z)}{\Gamma(\nu)
#' 2^{\nu - 1}}} with \eqn{z = 2\sqrt{\nu}\, h / \theta}, where \eqn{\theta} is
#' the range, \eqn{\nu} the smoothness and \eqn{K_\nu} the modified Bessel
#' function of the second kind. This is the parametrisation used by `geoR`.
#'
#' @param dmat A matrix of pairwise distances.
#' @param sigma2 Marginal variance (partial sill), positive.
#' @param theta Range parameter, positive.
#' @param nu Smoothness parameter, positive.
#' @return A covariance matrix of the same dimension as `dmat`.
#' @keywords internal
#' @noRd
matern_cov_matrix <- function(dmat, sigma2, theta, nu) {
  sigma2 <- max(as.numeric(sigma2), 1e-8)
  theta <- max(as.numeric(theta), 1e-8)
  nu <- max(as.numeric(nu), 1e-4)

  z <- 2 * sqrt(nu) * dmat / theta

  C <- matrix(0, nrow = nrow(dmat), ncol = ncol(dmat))
  positive <- z > 0
  C[positive] <- sigma2 * (z[positive]^nu) * besselK(z[positive], nu) /
    (gamma(nu) * 2^(nu - 1))
  diag(C) <- sigma2
  C[!is.finite(C)] <- 0
  C
}

#' Gaussian (restricted to constant mean) negative log-likelihood
#'
#' @param par Vector `c(log_sigma2, log_theta, log_nu, log_nugget)`.
#' @param r Residual vector.
#' @param dmat Distance matrix.
#' @return Negative log-likelihood (profiling out the constant mean), or a large
#'   finite penalty on numerical failure.
#' @keywords internal
#' @noRd
matern_negloglik <- function(par, r, dmat, nu = NULL) {
  sigma2 <- exp(par[1])
  theta <- exp(par[2])
  if (is.null(nu)) {
    nu <- exp(par[3])
    nugget <- exp(par[4])
  } else {
    nugget <- exp(par[3])
  }

  n <- length(r)
  C <- matern_cov_matrix(dmat, sigma2, theta, nu) + nugget * diag(n)

  ch <- tryCatch(chol(C), error = function(e) NULL)
  if (is.null(ch)) return(1e10)

  ## profile out the constant mean via generalised least squares
  ones <- rep(1, n)
  Cinv_one <- backsolve(ch, forwardsolve(t(ch), ones))
  Cinv_r <- backsolve(ch, forwardsolve(t(ch), r))
  mu <- sum(ones * Cinv_r) / sum(ones * Cinv_one)

  resid <- r - mu
  Cinv_resid <- backsolve(ch, forwardsolve(t(ch), resid))
  quad <- sum(resid * Cinv_resid)
  logdet <- 2 * sum(log(diag(ch)))

  val <- 0.5 * (logdet + quad + n * log(2 * pi))
  if (!is.finite(val)) 1e10 else val
}

#' Estimate Matern covariance parameters from spatial residuals
#'
#' Fits a Matern covariance with nugget to a residual field by maximum
#' likelihood. The default engine (`method = "mle"`) uses a self-contained
#' `stats::optim()` optimiser and has no external dependency; `method = "geoR"`
#' reproduces the reference analysis via [geoR::likfit()] when that package is
#' installed. On any numerical failure a moment-based fallback is returned so
#' that downstream estimation can proceed.
#'
#' By default the smoothness is **fixed** at `nu = 0.5` (an exponential
#' covariance). The free four-parameter Matern MLE is poorly identified on weak
#' residual fields and can send the smoothness estimate to a numerical boundary
#' (observed `nu` in the hundreds); fixing `nu` removes this instability and
#' matches the common exponential-field assumption. Pass `nu = NULL` to estimate
#' the smoothness freely (legacy behaviour).
#'
#' @param resid Numeric residual vector.
#' @param coords Two-column matrix of coordinates.
#' @param method Estimation engine, `"mle"` (default) or `"geoR"`.
#' @param nu Fixed Matern smoothness (default `0.5`, exponential); `NULL` to
#'   estimate it. Ignored by the `"geoR"` engine, which always estimates it.
#' @return A list with `sigma2`, `theta`, `nu` and `sigma2_eps` (nugget).
#' @keywords internal
#' @noRd
estimate_matern_params <- function(resid, coords, method = c("mle", "geoR"),
                                   nu = 0.5) {
  method <- match.arg(method)
  resid <- as.numeric(resid)

  vres <- stats::var(resid, na.rm = TRUE)
  if (!is.finite(vres) || vres <= 0) vres <- 1

  fallback <- list(sigma2 = 0.7 * vres, theta = 0.2,
                   nu = if (is.null(nu)) 0.5 else nu, sigma2_eps = 0.3 * vres)

  if (method == "geoR") {
    if (!requireNamespace("geoR", quietly = TRUE)) {
      warning("Package 'geoR' is not installed; falling back to method = 'mle'.",
              call. = FALSE)
      method <- "mle"
    } else {
      fit <- tryCatch({
        geodat <- geoR::as.geodata(
          data.frame(x = coords[, 1], y = coords[, 2], r = resid),
          coords.col = 1:2, data.col = 3)
        geoR::likfit(geodat, trend = "cte", cov.model = "matern",
                     ini.cov.pars = c(0.7 * vres, 0.2), nugget = 0.3 * vres,
                     kappa = 0.5, fix.nugget = FALSE, fix.kappa = FALSE,
                     messages = FALSE)
      }, error = function(e) NULL, warning = function(w) NULL)

      if (is.null(fit)) return(fallback)
      return(list(sigma2 = fit$cov.pars[1], theta = fit$cov.pars[2],
                  nu = fit$kappa, sigma2_eps = fit$nugget))
    }
  }

  dmat <- as.matrix(stats::dist(coords))

  if (is.null(nu)) {
    ## legacy: estimate smoothness too (four parameters)
    start <- log(c(0.7 * vres, 0.2, 0.5, 0.3 * vres))
    fit <- tryCatch(
      stats::optim(start, matern_negloglik, r = resid, dmat = dmat, nu = NULL,
                   method = "Nelder-Mead",
                   control = list(maxit = 500, reltol = 1e-8)),
      error = function(e) NULL)
    if (is.null(fit) || fit$convergence > 1) return(fallback)
    pars <- exp(fit$par)
    out <- list(sigma2 = pars[[1]], theta = pars[[2]], nu = pars[[3]],
                sigma2_eps = pars[[4]])
  } else {
    ## default: smoothness fixed at nu (three parameters, stable)
    start <- log(c(0.7 * vres, 0.2, 0.3 * vres))
    fit <- tryCatch(
      stats::optim(start, matern_negloglik, r = resid, dmat = dmat, nu = nu,
                   method = "Nelder-Mead",
                   control = list(maxit = 500, reltol = 1e-8)),
      error = function(e) NULL)
    if (is.null(fit) || fit$convergence > 1) return(fallback)
    pars <- exp(fit$par)
    out <- list(sigma2 = pars[[1]], theta = pars[[2]], nu = nu,
                sigma2_eps = pars[[3]])
  }
  if (!all(vapply(out, is.finite, logical(1)))) return(fallback)
  out
}
