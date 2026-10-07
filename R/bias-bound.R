#' Post-matching bias bound for an iDAPS fit
#'
#' Estimates an upper bound on the confounding bias of an [idaps()] matching
#' estimate, following the variogram bias bound. A matched-pair difference is a
#' spatial-differencing operation, so the estimator's bias decomposes into three
#' pieces controlled by the matched sample:
#' \deqn{|\mathrm{bias}| \le \|\theta_2\|\,\overline{\|\Delta X\|}
#'   + |\theta_3|\,\overline{|\Delta E|}
#'   + c\,|\theta_U|\,\overline{\sqrt{2\gamma_U(d)}},}
#' where \eqn{\Delta X, \Delta E} are within-pair covariate and neighbourhood-
#' exposure differences, \eqn{\gamma_U} is the (semi)variogram of the unmeasured
#' spatial confounder evaluated at the matched distances \eqn{d}, and \eqn{c} is
#' a selection constant (see `c_overlap`). The covariate and exposure
#' coefficients are taken from a working outcome model, and the confounding term
#' from the fitted spatial variogram of that model's residuals (whose sill
#' carries \eqn{\theta_U^2\,\mathrm{Var}(U)}). The confounding term shrinks as
#' the matched distances shrink, so a smaller caliper yields a smaller bound.
#'
#' @param object An `idaps_fit` returned by [idaps()] (supplies the matched
#'   pairs and neighbourhood exposure).
#' @param Y,Z,X,coords The same outcome, treatment, covariates and coordinates
#'   passed to [idaps()].
#' @param tau Exposure-kernel bandwidth, used only if `object` lacks a stored
#'   exposure (default `0.1`).
#' @param c_overlap Selection inflation constant on the confounding term
#'   (default `1.2`). Matching selects units with systematically larger
#'   confounder values, so the marginal variogram slightly understates the
#'   within-pair confounder differences; `c_overlap` \eqn{\approx 1.2} restores a
#'   finite-sample-valid bound (set to `1` for the population-level bound).
#' @param matern_method,matern_nu Passed to the residual variogram fit; see
#'   [recoverU()].
#'
#' @return An object of class `spaci_bias_bound`: a list with the covariate,
#'   exposure and confounding bound terms (`term_X`, `term_E`, `term_U`), their
#'   `total`, the mean matched distance, and the fitted variogram parameters.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 1)
#' fit <- idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, seed = 1)
#' bias_bound(fit, sim$Y, sim$Z, sim$X, sim$coords)
#' @export
bias_bound <- function(object, Y, Z, X, coords, tau = 0.1, c_overlap = 1.2,
                       matern_method = c("mle", "geoR"), matern_nu = 0.5) {
  matern_method <- match.arg(matern_method)
  if (!inherits(object, "idaps_fit")) {
    stop("`object` must be an idaps_fit returned by idaps().", call. = FALSE)
  }
  pairs <- object$extras$pairs
  if (is.null(pairs) || nrow(pairs) < 2) {
    stop("`object` has no matched pairs to bound.", call. = FALSE)
  }
  check_scalars(tau = tau)
  if (!is.numeric(c_overlap) || length(c_overlap) != 1L || c_overlap <= 0) {
    stop("`c_overlap` must be a single positive number.", call. = FALSE)
  }
  d <- validate_inputs(Y, Z, X, coords, require_X = TRUE, require_coords = TRUE)
  Y <- d$Y; Z <- d$Z; X <- d$X; coords <- d$coords
  xn <- colnames(X)

  E <- object$extras$E
  if (is.null(E)) E <- neighbourhood_exposure(coords, Z, tau = tau)$E

  ## working outcome model for the covariate / exposure coefficients + residuals
  dat <- data.frame(Y = Y, Z = Z, as.data.frame(X), G = E)
  fit <- stats::lm(stats::as.formula(
    paste("Y ~ Z +", paste(c(xn, "G"), collapse = " + "))), data = dat)
  co <- stats::coef(fit)
  th2 <- co[xn]; th3 <- unname(co["G"])
  ph <- estimate_matern_params(stats::residuals(fit), coords,
                               method = matern_method, nu = matern_nu)

  ti <- pairs[, 1]; ci <- pairs[, 2]
  dX <- X[ti, , drop = FALSE] - X[ci, , drop = FALSE]
  dE <- E[ti] - E[ci]
  d_pair <- sqrt(rowSums((coords[ti, , drop = FALSE] -
                          coords[ci, , drop = FALSE])^2))

  ## Matern correlation at the matched distances -> spatial semivariance
  z <- 2 * sqrt(ph$nu) * d_pair / ph$theta
  rho <- ifelse(z > 0,
                (z^ph$nu) * besselK(z, ph$nu) / (gamma(ph$nu) * 2^(ph$nu - 1)),
                1)
  rho[!is.finite(rho)] <- 0
  gU <- ph$sigma2 * (1 - rho)

  ## coordinatewise (l1) aggregation: each |theta_k| E|dX_k| is in outcome
  ## units, so the term is invariant to rescaling individual covariates (the
  ## Cauchy-Schwarz aggregate is not, and can be vacuous when covariate
  ## scales are heterogeneous)
  term_X <- sum(abs(th2) * colMeans(abs(dX)))
  term_E <- abs(th3) * mean(abs(dE))
  term_U <- c_overlap * mean(sqrt(2 * pmax(gU, 0)))

  structure(
    list(term_X = term_X, term_E = term_E, term_U = term_U,
         total = term_X + term_E + term_U, c_overlap = c_overlap,
         n_pairs = nrow(pairs), mean_distance = mean(d_pair),
         variogram = ph),
    class = "spaci_bias_bound")
}

#' Print method for a bias bound
#'
#' @param x A `spaci_bias_bound` object.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.spaci_bias_bound <- function(x, digits = 4, ...) {
  cat("Estimated post-matching bias bound (iDAPS)\n\n")
  cat("  covariate term     :", format(x$term_X, digits = digits), "\n")
  cat("  exposure term      :", format(x$term_E, digits = digits), "\n")
  cat("  confounding term   :", format(x$term_U, digits = digits),
      sprintf("(c_overlap = %g)\n", x$c_overlap))
  cat("  --------------------------------\n")
  cat("  total bias bound   :", format(x$total, digits = digits), "\n\n")
  cat("  matched pairs      :", x$n_pairs, "\n")
  cat("  mean matched dist. :", format(x$mean_distance, digits = digits), "\n")
  invisible(x)
}
