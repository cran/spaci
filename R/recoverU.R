#' Doubly robust ATT with a recovered spatial confounder (core engine)
#'
#' Shared implementation of the recoverU and recoverU+ estimators. An initial
#' outcome model `Y ~ Z + X + G` is fitted; a Matern field is estimated on its
#' residuals and used to recover the treatment-uncorrelated part of the spatial
#' confounder \eqn{U_R(s)} by generalised least squares (Equation 2.9). The
#' recovered confounder augments a propensity-score model and a control-outcome
#' model, which are combined in the doubly robust ATT estimator of Moodie et al.
#' (2018) / Equation (2.12). recoverU+ additionally includes the
#' neighbourhood-exposure term `G` in the PS and control-outcome models.
#'
#' @param Y Numeric outcome vector.
#' @param Z Binary treatment vector.
#' @param X Numeric covariate matrix.
#' @param E Numeric neighbourhood-exposure vector (`G`).
#' @param coords Two-column coordinate matrix.
#' @param d_space_raw Raw Euclidean distance matrix.
#' @param include_E_in_PS Logical; include `E` in the PS and control-outcome
#'   models (recoverU+ when `TRUE`, recoverU when `FALSE`).
#' @param matern_method Matern estimation engine, `"mle"` or `"geoR"`.
#' @return A list with `att`, `se`, `Uhat` (the recovered confounder) and `ps`
#'   (fitted propensity scores from the augmented model).
#' @keywords internal
#' @noRd
recoverU_core <- function(Y, Z, X, E, coords, d_space_raw,
                          include_E_in_PS = FALSE,
                          matern_method = "mle", matern_nu = 0.5) {
  n_obs <- length(Y)
  n1 <- sum(Z == 1)
  fail <- list(att = NA_real_, se = NA_real_, Uhat = rep(NA_real_, n_obs),
              ps = rep(NA_real_, n_obs))
  if (n1 < 2 || sum(Z == 0) < 2) return(fail)

  Xdf <- as.data.frame(X)
  xnames <- colnames(X)
  dat <- data.frame(Y = Y, Z = Z, Xdf, G = E)

  ## --- initial outcome model (always includes G) ----------------------------
  init_form <- stats::as.formula(
    paste("Y ~ Z +", paste(c(xnames, "G"), collapse = " + ")))
  fit_initial <- tryCatch(stats::lm(init_form, data = dat),
                          error = function(e) NULL)
  if (is.null(fit_initial)) return(fail)
  resid_initial <- stats::residuals(fit_initial)

  ## --- recover the spatial confounder via GLS (Eq. 2.9) ---------------------
  par_hat <- estimate_matern_params(resid_initial, coords, method = matern_method,
                                    nu = matern_nu)
  Sigma_U <- matern_cov_matrix(d_space_raw, par_hat$sigma2, par_hat$theta,
                               par_hat$nu)
  sigma2_eps_hat <- max(par_hat$sigma2_eps, 1e-8)
  V <- Sigma_U + sigma2_eps_hat * diag(n_obs) + 1e-6 * diag(n_obs)

  W <- stats::model.matrix(init_form, data = dat)

  theta_gls <- tryCatch({
    Vinv_W <- solve(V, W)
    Vinv_Y <- solve(V, Y)
    solve(t(W) %*% Vinv_W, t(W) %*% Vinv_Y)
  }, error = function(e) NULL)
  if (is.null(theta_gls)) return(fail)

  resid_gls <- as.numeric(Y - W %*% theta_gls)
  Uhat <- tryCatch(as.vector(Sigma_U %*% solve(V, resid_gls)),
                   error = function(e) rep(NA_real_, n_obs))
  if (any(!is.finite(Uhat))) return(fail)

  dat$Uhat <- safe_scale(Uhat)

  ## --- augmented PS and control-outcome models ------------------------------
  rhs <- c(xnames, "Uhat")
  if (include_E_in_PS) rhs <- c(xnames, "G", "Uhat")

  ps_form <- stats::as.formula(paste("Z ~", paste(rhs, collapse = " + ")))
  m0_form <- stats::as.formula(paste("Y ~", paste(rhs, collapse = " + ")))

  ps_fit <- tryCatch(
    suppressWarnings(stats::glm(ps_form, family = stats::binomial(),
                                data = dat,
                                control = stats::glm.control(maxit = 100))),
    error = function(e) NULL)
  if (is.null(ps_fit)) return(fail)
  ehat <- clip_ps(stats::fitted(ps_fit))

  control_dat <- dat[dat$Z == 0, , drop = FALSE]
  m0_fit <- tryCatch(stats::lm(m0_form, data = control_dat),
                     error = function(e) NULL)
  if (is.null(m0_fit)) return(fail)

  m0hat <- tryCatch(as.numeric(stats::predict(m0_fit, newdata = dat)),
                    error = function(e) rep(NA_real_, n_obs))
  if (any(!is.finite(m0hat))) return(fail)

  ## --- doubly robust ATT (Eq. 2.12) -----------------------------------------
  w <- ehat / (1 - ehat)
  psi <- (Z - (1 - Z) * w) * (Y - m0hat)
  att_hat <- sum(psi) / n1

  pi_hat <- n1 / n_obs
  infl <- psi / pi_hat - att_hat
  se_hat <- stats::sd(infl, na.rm = TRUE) / sqrt(n_obs)

  list(att = att_hat, se = se_hat, Uhat = Uhat, infl = infl, ps = ehat)
}

#' recoverU: doubly robust ATT with a recovered spatial confounder
#'
#' Doubly robust ATT estimator that adjusts for spatial confounding through a
#' partially recovered spatial confounder (Pokal et al., 2023), but does not
#' adjust the propensity-score / control-outcome models for spatial
#' interference. Provided as a comparator to [recoverUplus()].
#'
#' @param Y Numeric outcome vector.
#' @param Z Binary treatment vector (0/1).
#' @param X Covariate matrix or data frame.
#' @param coords Two-column matrix or data frame of coordinates.
#' @param tau Positive spatial decay parameter for the exposure kernel used in
#'   the initial outcome model (default `0.1`).
#' @param normalize Logical; row-normalise the exposure kernel (default `TRUE`).
#' @param matern_method Matern estimation engine: `"mle"` (default, no external
#'   dependency) or `"geoR"` (reproduces the reference analysis).
#' @param matern_nu Fixed Matern smoothness for the `"mle"` engine (default
#'   `0.5`, an exponential covariance, which is numerically stable); pass `NULL`
#'   to estimate the smoothness freely. Ignored by the `"geoR"` engine.
#' @param level Confidence level for the reported interval.
#'
#' @return An `idaps_fit` object.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 1)
#' recoverU(sim$Y, sim$Z, sim$X, sim$coords)
#' @export
recoverU <- function(Y, Z, X, coords, tau = 0.1, normalize = TRUE,
                     matern_method = c("mle", "geoR"), matern_nu = 0.5,
                     level = 0.95) {
  matern_method <- match.arg(matern_method)
  check_scalars(tau = tau, level = level)
  d <- validate_inputs(Y, Z, X, coords, require_X = TRUE, require_coords = TRUE)
  ex <- neighbourhood_exposure(d$coords, d$Z, tau = tau, normalize = normalize)
  d_space_raw <- as.matrix(stats::dist(d$coords))

  res <- recoverU_core(d$Y, d$Z, d$X, ex$E, d$coords, d_space_raw,
                       include_E_in_PS = FALSE, matern_method = matern_method,
                       matern_nu = matern_nu)
  new_idaps_fit("recoverU", res$att, res$se,
                extras = list(E = ex$E, Uhat = res$Uhat, infl = res$infl,
                              ps = res$ps, coords = d$coords), level = level)
}

#' recoverU+: doubly robust ATT under spatial confounding and interference
#'
#' The recoverU+ estimator augments the doubly robust ATT with both a partially
#' recovered spatial confounder and a neighbourhood-exposure term, so that the
#' propensity-score and control-outcome models adjust for spatial confounding
#' and spatial interference simultaneously (Equations 2.10-2.12). Compared with
#' [recoverU()], the neighbourhood exposure `G` is included in the PS and
#' control-outcome models.
#'
#' @inheritParams recoverU
#'
#' @return An `idaps_fit` object.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 1)
#' recoverUplus(sim$Y, sim$Z, sim$X, sim$coords)
#' @export
recoverUplus <- function(Y, Z, X, coords, tau = 0.1, normalize = TRUE,
                         matern_method = c("mle", "geoR"), matern_nu = 0.5,
                         level = 0.95) {
  matern_method <- match.arg(matern_method)
  check_scalars(tau = tau, level = level)
  d <- validate_inputs(Y, Z, X, coords, require_X = TRUE, require_coords = TRUE)
  ex <- neighbourhood_exposure(d$coords, d$Z, tau = tau, normalize = normalize)
  d_space_raw <- as.matrix(stats::dist(d$coords))

  res <- recoverU_core(d$Y, d$Z, d$X, ex$E, d$coords, d_space_raw,
                       include_E_in_PS = TRUE, matern_method = matern_method,
                       matern_nu = matern_nu)
  new_idaps_fit("recoverU+", res$att, res$se,
                extras = list(E = ex$E, Uhat = res$Uhat, infl = res$infl,
                              ps = res$ps, coords = d$coords), level = level)
}
