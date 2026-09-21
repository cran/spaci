# idaps 0.1.0

* First release. Implements the `idaps()` (distance-adjusted propensity score
  with interference) and `recoverUplus()` (doubly robust with recovered spatial
  confounder and neighbourhood exposure) estimators, together with the
  `naive_ps()`, `daps()` and `recoverU()` comparators.
* `simulate_spatial_causal()` generates data with joint spatial confounding and
  interference; `spatial_ate()` runs all estimators and collects the results.
* The recovered-confounder step uses a self-contained Matérn maximum-likelihood
  fit by default (`matern_method = "mle"`), with an optional `geoR` engine for
  faithful reproduction of the report.
