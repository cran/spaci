#' Simulate spatial data with confounding and interference
#'
#' Draws one data set from the data-generating process of the report
#' (Equation 2.13 and Section 2.3): an unmeasured spatial confounder \eqn{U(s)}
#' drawn from an exponential Gaussian random field drives both treatment and
#' outcome, and the outcome further depends on a neighbourhood-exposure term so
#' that spatial interference is present. The true ATT is `true_att`.
#'
#' The Gaussian random field is generated directly from its covariance matrix
#' (exponential covariance, i.e. Matern with smoothness 1/2), so the simulator
#' has no external dependency.
#'
#' @param n Number of units.
#' @param true_att True average treatment effect on the treated.
#' @param beta0,beta1,beta2 Outcome-model intercept and covariate coefficients.
#' @param theta_spatial Coefficient of the spatial confounder in the outcome.
#' @param gamma_interference Coefficient of neighbourhood exposure in the
#'   outcome (interference strength).
#' @param sigma_eps Outcome noise standard deviation.
#' @param delta_u Strength of the spatial confounder in the treatment model
#'   (controls the degree of spatial confounding).
#' @param u_phi Range parameter of the confounder's exponential covariance.
#' @param tau_exp Spatial decay parameter of the exposure kernel.
#' @param normalize Logical; row-normalise the exposure kernel (default `TRUE`).
#' @param seed Optional integer seed.
#'
#' @return A list with `Y`, `Z`, `X` (matrix with columns `X1`, `X2`), `coords`,
#'   the exposure `E`, the latent confounder `U`, and `true_att`.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 42)
#' str(sim)
#' @export
simulate_spatial_causal <- function(n = 150, true_att = 2.0,
                                     beta0 = 2.5, beta1 = 1.0, beta2 = 0.5,
                                     theta_spatial = 0.4,
                                     gamma_interference = 1.5,
                                     sigma_eps = 1.0, delta_u = 2.0,
                                     u_phi = 0.2, tau_exp = 0.1,
                                     normalize = TRUE, seed = NULL) {
  if (!is.null(seed)) {
    old <- .save_seed(); on.exit(.restore_seed(old), add = TRUE); set.seed(seed)
  }

  coords <- cbind(stats::runif(n), stats::runif(n))
  d_space_raw <- as.matrix(stats::dist(coords))

  X1 <- stats::rnorm(n)
  X2 <- stats::rbinom(n, 1, 0.5)
  X <- cbind(X1 = X1, X2 = X2)

  ## exponential Gaussian random field for the spatial confounder
  Sigma <- exp(-d_space_raw / u_phi)
  L <- chol(Sigma + 1e-8 * diag(n))
  U <- as.vector(crossprod(L, stats::rnorm(n)))
  U <- as.vector(scale(U))
  U[!is.finite(U)] <- 0

  ps_true <- clip_ps(stats::plogis(0.1 + 0.1 * X1 + 0.2 * X2 + delta_u * U))
  Z <- stats::rbinom(n, 1, ps_true)
  if (length(unique(Z)) < 2) Z <- stats::rbinom(n, 1, 0.5)

  ex <- neighbourhood_exposure(coords, Z, tau = tau_exp, normalize = normalize)
  E <- ex$E

  eps <- stats::rnorm(n, 0, sigma_eps)
  Y <- beta0 + true_att * Z + beta1 * X1 + beta2 * X2 +
    theta_spatial * U + gamma_interference * E + eps

  list(Y = Y, Z = Z, X = X, coords = coords, E = E, U = U,
       true_att = true_att)
}

#' Estimate the ATT with all available methods
#'
#' Convenience wrapper that runs the naive propensity score, DAPS, iDAPS,
#' recoverU and recoverU+ estimators on the same data and collects the point
#' estimates, standard errors and confidence intervals in a single data frame.
#'
#' @inheritParams idaps
#' @param matern_method Matern estimation engine for the recoverU family,
#'   `"mle"` (default) or `"geoR"`.
#' @param methods Character vector selecting which estimators to run. Defaults
#'   to all five.
#'
#' @return A data frame with one row per method (`Method`, `ATT`, `SE`,
#'   `Lower`, `Upper`), with the fitted objects attached as the `"fits"`
#'   attribute.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 1)
#' spatial_ate(sim$Y, sim$Z, sim$X, sim$coords, seed = 1)
#' @export
spatial_ate <- function(Y, Z, X, coords, tau = 0.1, caliper = 0.25,
                        matern_method = c("mle", "geoR"), level = 0.95,
                        seed = NULL,
                        methods = c("Naive PS", "DAPS", "iDAPS",
                                    "recoverU", "recoverU+")) {
  matern_method <- match.arg(matern_method)
  methods <- match.arg(methods, several.ok = TRUE)

  runners <- list(
    "Naive PS"  = function() naive_ps(Y, Z, X, caliper = caliper,
                                      seed = seed, level = level),
    "DAPS"      = function() daps(Y, Z, X, coords, caliper = caliper,
                                  seed = seed, level = level),
    "iDAPS"     = function() idaps(Y, Z, X, coords, tau = tau,
                                   caliper = caliper, seed = seed,
                                   level = level),
    "recoverU"  = function() recoverU(Y, Z, X, coords, tau = tau,
                                      matern_method = matern_method,
                                      level = level),
    "recoverU+" = function() recoverUplus(Y, Z, X, coords, tau = tau,
                                          matern_method = matern_method,
                                          level = level)
  )

  fits <- lapply(methods, function(m) runners[[m]]())
  names(fits) <- methods

  out <- data.frame(
    Method = methods,
    ATT = vapply(fits, function(f) f$att, numeric(1)),
    SE = vapply(fits, function(f) f$se, numeric(1)),
    Lower = vapply(fits, function(f) f$ci[["lower"]], numeric(1)),
    Upper = vapply(fits, function(f) f$ci[["upper"]], numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
  attr(out, "fits") <- fits
  out
}
