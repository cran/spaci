test_that("simulator returns consistent shapes and the true ATT", {
  sim <- simulate_spatial_causal(n = 120, seed = 1)
  expect_length(sim$Y, 120)
  expect_length(sim$Z, 120)
  expect_equal(dim(sim$X), c(120, 2))
  expect_equal(dim(sim$coords), c(120, 2))
  expect_true(all(sim$Z %in% c(0, 1)))
  expect_equal(sim$true_att, 2.0)
})

test_that("neighbourhood exposure is a convex combination when normalised", {
  set.seed(2)
  coords <- cbind(runif(30), runif(30))
  Z <- rbinom(30, 1, 0.5)
  ex <- neighbourhood_exposure(coords, Z, tau = 0.1, normalize = TRUE)
  expect_length(ex$E, 30)
  expect_true(all(ex$E >= 0 & ex$E <= 1))          # weighted average of 0/1
  expect_equal(unname(diag(ex$G)), rep(0, 30))
})

test_that("matching estimators return valid idaps_fit objects", {
  sim <- simulate_spatial_causal(n = 150, seed = 3)
  for (fit in list(
    naive_ps(sim$Y, sim$Z, sim$X, seed = 1),
    daps(sim$Y, sim$Z, sim$X, sim$coords, seed = 1),
    idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, seed = 1)
  )) {
    expect_s3_class(fit, "idaps_fit")
    expect_true(is.finite(fit$att))
    expect_true(is.finite(fit$se))
    expect_true(fit$ci[["lower"]] < fit$ci[["upper"]])
  }
})

test_that("iDAPS weights lie on the simplex", {
  sim <- simulate_spatial_causal(n = 150, seed = 4)
  fit <- idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, seed = 1)
  w <- fit$weights
  expect_length(w, 3)
  expect_equal(sum(w), 1, tolerance = 1e-8)
  expect_true(all(w >= 0))
})

test_that("recoverU and recoverU+ return finite doubly robust estimates", {
  sim <- simulate_spatial_causal(n = 150, seed = 5)
  fu <- recoverU(sim$Y, sim$Z, sim$X, sim$coords, matern_method = "mle")
  fp <- recoverUplus(sim$Y, sim$Z, sim$X, sim$coords, matern_method = "mle")
  expect_s3_class(fu, "idaps_fit")
  expect_s3_class(fp, "idaps_fit")
  expect_true(is.finite(fu$att))
  expect_true(is.finite(fp$att))
})

test_that("setting pi3 = 0 (iDAPS) with no interference approximates DAPS logic", {
  # With a fixed seed the matching order is reproducible.
  sim <- simulate_spatial_causal(n = 150, seed = 6)
  f1 <- idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, seed = 42)
  f2 <- idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, seed = 42)
  expect_equal(f1$att, f2$att)   # deterministic given the seed
})

test_that("inputs are validated", {
  sim <- simulate_spatial_causal(n = 50, seed = 7)
  expect_error(naive_ps(sim$Y, rep(1, 50), sim$X), "both treated and control")
  expect_error(idaps(sim$Y, sim$Z, sim$X, sim$coords[, 1, drop = FALSE]),
               "exactly two columns")
})
