test_that("scalar arguments are validated with clear messages", {
  sim <- simulate_spatial_causal(n = 80, seed = 1)
  expect_error(naive_ps(sim$Y, sim$Z, sim$X, caliper = -1), "caliper")
  expect_error(naive_ps(sim$Y, sim$Z, sim$X, caliper = 0), "caliper")
  expect_error(idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0), "tau")
  expect_error(idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = -0.5), "tau")
  expect_error(recoverU(sim$Y, sim$Z, sim$X, sim$coords, level = 1.2), "level")
  expect_error(idaps(sim$Y, sim$Z, sim$X, sim$coords, pi_grid_step = 2),
               "grid step")
})

test_that("missing / malformed data are rejected", {
  sim <- simulate_spatial_causal(n = 80, seed = 2)
  y <- sim$Y; y[1] <- NA
  expect_error(naive_ps(y, sim$Z, sim$X), "missing")
  expect_error(daps(sim$Y, sim$Z, NULL, sim$coords), "must be supplied")
  expect_error(idaps(sim$Y, sim$Z, sim$X, NULL), "must be supplied")
  expect_error(idaps(sim$Y, sim$Z, sim$X, sim$coords[, 1, drop = FALSE]),
               "two columns")
  x <- sim$X; x[1, 1] <- Inf
  expect_error(naive_ps(sim$Y, sim$Z, x), "finite")
})

test_that("degenerate treatment is rejected", {
  sim <- simulate_spatial_causal(n = 80, seed = 3)
  expect_error(naive_ps(sim$Y, rep(1, 80), sim$X), "both treated and control")
})

test_that("plot_ate accepts spatial_ate output and a list of fits", {
  sim <- simulate_spatial_causal(n = 120, seed = 4)
  res <- spatial_ate(sim$Y, sim$Z, sim$X, sim$coords, seed = 1)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp); on.exit(unlink(tmp))
  out <- plot_ate(res, true_att = sim$true_att)
  fit <- idaps(sim$Y, sim$Z, sim$X, sim$coords, seed = 1)
  plot(fit)                                # single-fit plot method
  grDevices::dev.off()
  expect_true(all(c("Method", "ATT", "Lower", "Upper") %in% names(out)))
  expect_equal(nrow(out), 5)
})

test_that("as_ate_df rejects unsupported input", {
  expect_error(plot_ate(data.frame(a = 1)), "must have columns")
  expect_error(plot_ate(42), "spatial_ate")
})
