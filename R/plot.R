#' Forest plot of ATT estimates across methods
#'
#' Draws a forest (caterpillar) plot of the estimated ATT and its confidence
#' interval for each method, with a reference line at zero and, optionally, at
#' the true effect. This reproduces the style of Figure 2.5 of the report and is
#' the recommended way to compare methods visually.
#'
#' @param x Either the data frame returned by [spatial_ate()] (columns
#'   `Method`, `ATT`, `Lower`, `Upper`) or a named list of
#'   `idaps_fit` objects.
#' @param true_att Optional true effect; drawn as a dashed vertical reference
#'   line when supplied.
#' @param null_line Position of the "no effect" reference line (default `0`);
#'   set to `NA` to omit.
#' @param col,pt_col Colours for the confidence intervals and point estimates.
#' @param xlab,main Axis label and title.
#' @param ... Passed to [graphics::plot()].
#'
#' @return The plotted data frame, invisibly.
#' @examples
#' sim <- simulate_spatial_causal(n = 50, seed = 1)
#' res <- spatial_ate(sim$Y, sim$Z, sim$X, sim$coords, seed = 1)
#' plot_ate(res, true_att = sim$true_att)
#' @export
plot_ate <- function(x, true_att = NULL, null_line = 0,
                     col = "#4C72B0", pt_col = "#1A1A1A",
                     xlab = "Average treatment effect on the treated (ATT)",
                     main = "Estimated ATT and confidence intervals", ...) {
  df <- as_ate_df(x)
  if (nrow(df) == 0L) stop("No estimates to plot.", call. = FALSE)

  df <- df[rev(seq_len(nrow(df))), , drop = FALSE]  # first method at the top
  yy <- seq_len(nrow(df))

  finite_vals <- c(df$Lower, df$Upper, df$ATT, true_att, null_line)
  finite_vals <- finite_vals[is.finite(finite_vals)]
  xlim <- range(finite_vals)
  pad <- 0.05 * diff(xlim)
  if (pad == 0) pad <- 0.5
  xlim <- xlim + c(-pad, pad)

  op <- graphics::par(mar = c(4.5, 7, 3, 1))
  on.exit(graphics::par(op), add = TRUE)

  graphics::plot(df$ATT, yy, xlim = xlim, ylim = c(0.5, nrow(df) + 0.5),
                 pch = NA, yaxt = "n", ylab = "", xlab = xlab, main = main,
                 bty = "n", ...)
  graphics::axis(2, at = yy, labels = df$Method, las = 1, tick = FALSE)

  if (!is.na(null_line)) {
    graphics::abline(v = null_line, lty = 2, col = "grey50")
  }
  if (!is.null(true_att) && is.finite(true_att)) {
    graphics::abline(v = true_att, lty = 3, col = "#C44E52", lwd = 2)
    graphics::mtext(sprintf("true = %.3g", true_att), side = 3, adj = 1,
                    col = "#C44E52", cex = 0.8)
  }

  graphics::segments(df$Lower, yy, df$Upper, yy, col = col, lwd = 2)
  ends <- 0.12
  graphics::segments(df$Lower, yy - ends, df$Lower, yy + ends, col = col, lwd = 2)
  graphics::segments(df$Upper, yy - ends, df$Upper, yy + ends, col = col, lwd = 2)
  graphics::points(df$ATT, yy, pch = 19, col = pt_col, cex = 1.2)
  graphics::text(df$ATT, yy + 0.28, labels = formatC(df$ATT, format = "f",
                 digits = 2), cex = 0.8, col = pt_col)

  invisible(df[rev(seq_len(nrow(df))), , drop = FALSE])
}

#' Coerce supported inputs to a plotting data frame
#'
#' @param x A [spatial_ate()] data frame or a list of `idaps_fit` objects.
#' @return A data frame with `Method`, `ATT`, `Lower`, `Upper`.
#' @keywords internal
#' @noRd
as_ate_df <- function(x) {
  if (is.data.frame(x)) {
    need <- c("Method", "ATT", "Lower", "Upper")
    if (!all(need %in% names(x))) {
      stop("Data frame must have columns: ", paste(need, collapse = ", "),
           ".", call. = FALSE)
    }
    return(x[, need, drop = FALSE])
  }
  if (inherits(x, "idaps_fit")) x <- list(x)
  if (is.list(x) && all(vapply(x, inherits, logical(1), "idaps_fit"))) {
    return(data.frame(
      Method = vapply(x, function(f) f$method, character(1)),
      ATT = vapply(x, function(f) f$att, numeric(1)),
      Lower = vapply(x, function(f) f$ci[["lower"]], numeric(1)),
      Upper = vapply(x, function(f) f$ci[["upper"]], numeric(1)),
      stringsAsFactors = FALSE
    ))
  }
  stop("`x` must be a spatial_ate() data frame or a list of idaps_fit objects.",
       call. = FALSE)
}

#' Plot a single estimate
#'
#' @param x An `idaps_fit` object.
#' @param ... Passed to [plot_ate()].
#' @return The plotted data frame, invisibly.
#' @export
plot.idaps_fit <- function(x, ...) {
  plot_ate(list(x), ...)
}
