#' Construct an idaps fit object
#'
#' @param method Character label for the estimator.
#' @param att Point estimate of the ATT.
#' @param se Standard error.
#' @param n_match,n_drop Matched / dropped treated counts (matching estimators).
#' @param weights Named vector of tuning weights (matching estimators).
#' @param extras List of additional diagnostics.
#' @param level Confidence level for the interval.
#' @return An object of class `idaps_fit`.
#' @keywords internal
#' @noRd
new_idaps_fit <- function(method, att, se, n_match = NA_integer_,
                          n_drop = NA_integer_, weights = NULL,
                          extras = list(), level = 0.95) {
  z <- stats::qnorm(1 - (1 - level) / 2)
  ci <- c(lower = att - z * se, upper = att + z * se)
  structure(
    list(method = method, estimand = "ATT", att = att, se = se,
         ci = ci, level = level, n_match = n_match, n_drop = n_drop,
         weights = weights, extras = extras),
    class = "idaps_fit"
  )
}

#' Print method for idaps fit objects
#'
#' @param x An `idaps_fit` object.
#' @param digits Number of significant digits.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.idaps_fit <- function(x, digits = 4, ...) {
  cat("Spatial causal effect estimate (", x$method, ")\n", sep = "")
  cat("Estimand: average treatment effect on the treated (ATT)\n\n")
  cat("  ATT   =", format(x$att, digits = digits), "\n")
  cat("  SE    =", format(x$se, digits = digits), "\n")
  cat(sprintf("  %g%% CI = [%s, %s]\n", 100 * x$level,
              format(x$ci[["lower"]], digits = digits),
              format(x$ci[["upper"]], digits = digits)))
  if (!is.null(x$weights)) {
    cat("\n  Tuning weights:",
        paste(sprintf("%s=%s", names(x$weights),
                      format(x$weights, digits = digits)), collapse = ", "),
        "\n")
  }
  if (is.finite(x$n_match)) {
    cat("  Matched treated units:", x$n_match,
        "(dropped:", x$n_drop, ")\n")
  }
  invisible(x)
}

#' Fit the naive propensity-score model
#'
#' @param Z Binary treatment vector.
#' @param X Numeric covariate matrix.
#' @return Clipped propensity scores.
#' @keywords internal
#' @noRd
estimate_ps <- function(Z, X) {
  dat <- data.frame(Z = Z, as.data.frame(X))
  form <- stats::as.formula(
    paste("Z ~", paste(colnames(X), collapse = " + ")))
  fit <- suppressWarnings(
    stats::glm(form, family = stats::binomial(link = "logit"), data = dat))
  clip_ps(stats::fitted(fit))
}
