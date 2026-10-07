# spaci 0.2.0

* Spatial-dependence-valid inference: `vcov_hac()` (Conley spatial HAC variance
  for `recoverU()`/`recoverUplus()` fits, with a data-driven bandwidth),
  `boot_spatial()` (spatial block bootstrap that re-runs the whole pipeline on
  resampled blocks), and `rand_test()` (conditional randomization test for the
  sharp null of no direct effect). The default i.i.d. standard error is
  anticonservative under spatial dependence — these give wider, better-calibrated
  intervals. `recoverU`/`recoverU+` fits now carry their influence values.
* The `"mle"` Matérn recovery engine now **fixes the smoothness at `nu = 0.5`
  (exponential) by default**, controlled by the new `matern_nu` argument of
  `recoverU()` / `recoverUplus()`. The free four-parameter fit is poorly
  identified on weak residual fields and could send the smoothness estimate to a
  numerical boundary; pass `matern_nu = NULL` for the previous behaviour. The
  `"geoR"` engine is unchanged.
* New `bias_bound()`: an estimated post-matching bias bound for an `idaps()` fit
  (covariate, exposure and variogram-based confounding terms), with a
  `print()` method.
* Matching estimators (`idaps()`, `daps()`, `naive_ps()`) gain a
  `match_method` argument: `"greedy"` (default, as before) or `"optimal"`, a
  deterministic 1:1 assignment via `clue::solve_LSAP()` that needs no seed.
* `boot_spatial()` and `rand_test()` restore the caller's random-number state
  when `seed` is supplied, as the other estimators already do.

# spaci 0.1.1

* First CRAN release.
* `simulate_spatial_causal()` now defaults to `n = 150` (previously 250).
* Examples run on smaller simulated data sets to keep check times short.

# spaci 0.1.0

* First release. Implements the `idaps()` (distance-adjusted propensity score
  with interference) and `recoverUplus()` (doubly robust with recovered spatial
  confounder and neighbourhood exposure) estimators, together with the
  `naive_ps()`, `daps()` and `recoverU()` comparators.
* `simulate_spatial_causal()` generates data with joint spatial confounding and
  interference; `spatial_ate()` runs all estimators and collects the results.
* The recovered-confounder step uses a self-contained Matérn maximum-likelihood
  fit by default (`matern_method = "mle"`), with an optional `geoR` engine for
  faithful reproduction of the report.
