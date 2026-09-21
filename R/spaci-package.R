#' spaci: Causal effect estimation under spatial confounding and interference
#'
#' The spaci package implements two unified methods for estimating the average
#' treatment effect on the treated (ATT) from spatial observational data when
#' spatial confounding and spatial interference are present simultaneously:
#'
#' \describe{
#'   \item{[idaps()]}{Distance-adjusted propensity score with interference: a
#'     matching estimator using a data-driven composite of propensity-score
#'     distance, spatial proximity and neighbourhood-exposure distance.}
#'   \item{[recoverUplus()]}{A doubly robust estimator that augments the
#'     propensity-score and control-outcome models with a partially recovered
#'     spatial confounder and a neighbourhood-exposure term.}
#' }
#'
#' The naive propensity score ([naive_ps()]), DAPS ([daps()]) and recoverU
#' ([recoverU()]) comparators, a data simulator ([simulate_spatial_causal()])
#' and an all-methods wrapper ([spatial_ate()]) are also provided.
#'
#' @keywords internal
"_PACKAGE"
