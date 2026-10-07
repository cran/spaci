test_that("vcov_hac returns a valid HAC variance for DR fits", {
  sim <- simulate_spatial_causal(n = 150, delta_u = 2, theta_spatial = 1.5, seed = 1)
  fit <- recoverUplus(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1)
  hac <- vcov_hac(fit)
  expect_s3_class(hac, "spaci_vcov")
  expect_true(is.finite(hac$se) && hac$se >= 0)
  expect_true(hac$variance >= 0)
  expect_true(hac$ci[["lower"]] < hac$ci[["upper"]])
  expect_true(is.finite(hac$bandwidth) && hac$bandwidth > 0)
  # a fixed numeric bandwidth also works
  expect_s3_class(vcov_hac(fit, bandwidth = 0.2, kernel = "uniform"), "spaci_vcov")
})

test_that("vcov_hac rejects matching fits (no influence values)", {
  sim <- simulate_spatial_causal(n = 120, seed = 2)
  fit <- idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, seed = 1)
  expect_error(vcov_hac(fit), "influence values")
})

test_that("boot_spatial returns a percentile interval", {
  sim <- simulate_spatial_causal(n = 120, delta_u = 2, theta_spatial = 1.5, seed = 3)
  bb <- boot_spatial(sim$Y, sim$Z, sim$X, sim$coords, method = "recoverUplus",
                     B = 25, seed = 1)
  expect_s3_class(bb, "spaci_boot")
  expect_true(is.finite(bb$att))
  expect_true(bb$ci[["lower"]] <= bb$ci[["upper"]])
  expect_gt(bb$B, 0)
})

test_that("boot_spatial works for a matching estimator too", {
  sim <- simulate_spatial_causal(n = 120, seed = 4)
  bb <- boot_spatial(sim$Y, sim$Z, sim$X, sim$coords, method = "idaps",
                     B = 20, seed = 1)
  expect_s3_class(bb, "spaci_boot")
})

test_that("rand_test returns a p-value in [0,1]", {
  sim <- simulate_spatial_causal(n = 120, delta_u = 2, theta_spatial = 1.5, seed = 5)
  rt <- rand_test(sim$Y, sim$Z, sim$X, sim$coords, method = "recoverUplus",
                  R = 25, seed = 1)
  expect_s3_class(rt, "spaci_randtest")
  expect_true(rt$p_value >= 0 && rt$p_value <= 1)
  expect_true(is.finite(rt$statistic))
})
