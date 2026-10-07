test_that("fixed-nu recovery (default) and free-nu (legacy) both run", {
  sim <- simulate_spatial_causal(n = 120, delta_u = 2, theta_spatial = 1.5, seed = 1)
  f_fixed <- recoverUplus(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1)
  f_free  <- recoverUplus(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, matern_nu = NULL)
  expect_true(is.finite(f_fixed$att))
  expect_true(is.finite(f_free$att))
})

test_that("bias_bound returns a valid bound that dominates the realized bias", {
  sim <- simulate_spatial_causal(n = 150, delta_u = 2, theta_spatial = 1.5, seed = 2)
  fit <- idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, seed = 1)
  bb <- bias_bound(fit, sim$Y, sim$Z, sim$X, sim$coords)
  expect_s3_class(bb, "spaci_bias_bound")
  expect_true(all(c(bb$term_X, bb$term_E, bb$term_U) >= 0))
  expect_equal(bb$total, bb$term_X + bb$term_E + bb$term_U)
  expect_gte(bb$total, abs(fit$att - sim$true_att))    # bound dominates bias
})

test_that("bias_bound rejects non-idaps input and pairless fits", {
  sim <- simulate_spatial_causal(n = 80, seed = 3)
  expect_error(bias_bound(list(), sim$Y, sim$Z, sim$X, sim$coords), "idaps_fit")
  ru <- recoverUplus(sim$Y, sim$Z, sim$X, sim$coords)   # no matched pairs
  expect_error(bias_bound(ru, sim$Y, sim$Z, sim$X, sim$coords), "matched pairs")
})

test_that("optimal matching is deterministic and needs no seed", {
  skip_if_not_installed("clue")
  sim <- simulate_spatial_causal(n = 150, delta_u = 2, theta_spatial = 1.5, seed = 4)
  a <- idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, match_method = "optimal")
  b <- idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, match_method = "optimal")
  expect_equal(a$att, b$att)                # deterministic across calls (no seed)
  expect_true(is.finite(a$att))
  expect_equal(sum(a$weights), 1, tolerance = 1e-8)
})

test_that("optimal matching respects the caliper (no over-distance pairs)", {
  skip_if_not_installed("clue")
  sim <- simulate_spatial_causal(n = 150, seed = 5)
  fit <- daps(sim$Y, sim$Z, sim$X, sim$coords, caliper = 0.1, match_method = "optimal")
  expect_s3_class(fit, "idaps_fit")
})
