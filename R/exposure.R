#' Neighbourhood exposure via a spatial kernel
#'
#' Constructs the neighbourhood-exposure mapping
#' \eqn{E_i = \sum_{j \ne i} G_{ij} A_j} that summarises spatial interference,
#' where \eqn{G_{ij}} is an exponential kernel of the distance between units
#' \eqn{i} and \eqn{j}, \eqn{G_{ij} = \exp(-d_{ij} / \tau)} with
#' \eqn{G_{ii} = 0}. This is the exposure used by [idaps()] and the
#' recoverU family.
#'
#' @param coords A two-column matrix or data frame of spatial coordinates.
#' @param Z Binary treatment vector (0/1), one entry per row of `coords`.
#' @param tau Positive spatial decay (bandwidth) parameter of the exponential
#'   kernel. Smaller values concentrate exposure on nearer neighbours.
#' @param normalize Logical; if `TRUE` (default) each row of the kernel is
#'   normalised to sum to one, so that exposure is a weighted average of
#'   neighbours' treatment. If `FALSE`, the raw kernel weights are used, matching
#'   the unnormalised definition in Equation (2.4) of the report.
#'
#' @return A list with components
#'   \describe{
#'     \item{E}{Numeric vector of neighbourhood exposures, one per unit.}
#'     \item{G}{The (possibly row-normalised) kernel matrix with zero diagonal.}
#'   }
#' @examples
#' set.seed(1)
#' coords <- cbind(runif(20), runif(20))
#' Z <- rbinom(20, 1, 0.5)
#' ex <- neighbourhood_exposure(coords, Z, tau = 0.1)
#' head(ex$E)
#' @export
neighbourhood_exposure <- function(coords, Z, tau = 0.1, normalize = TRUE) {
  coords <- as.matrix(coords)
  Z <- as.numeric(Z)
  if (nrow(coords) != length(Z)) {
    stop("`coords` and `Z` must describe the same number of units.",
         call. = FALSE)
  }
  if (!is.numeric(tau) || length(tau) != 1L || tau <= 0) {
    stop("`tau` must be a single positive number.", call. = FALSE)
  }

  d_raw <- as.matrix(stats::dist(coords))
  G <- exp(-d_raw / tau)
  diag(G) <- 0

  if (normalize) {
    row_sums <- rowSums(G)
    row_sums[row_sums == 0] <- 1
    G <- G / row_sums
  }

  E <- as.vector(G %*% Z)
  list(E = E, G = G)
}

#' Pairwise distance components for the composite matching metric
#'
#' Builds the three normalised distance matrices that make up the iDAPS
#' composite metric of Equation (2.3): the propensity-score distance
#' \eqn{D^{PS}_{ij} = |p(X_i) - p(X_j)|}, the spatial (Euclidean) distance
#' \eqn{D^{Spatial}_{ij}}, and the neighbourhood-exposure distance
#' \eqn{D^{Interference}_{ij} = |E_i - E_j|}. Each is rescaled to `[0, 1]` so
#' that no single component dominates the composite.
#'
#' @param ps Numeric vector of estimated propensity scores.
#' @param coords Two-column matrix of spatial coordinates.
#' @param E Numeric vector of neighbourhood exposures (see
#'   [neighbourhood_exposure()]).
#'
#' @return A list with the normalised matrices `D_ps`, `d_space`, `d_exp`
#'   (each with zero diagonal) and the raw Euclidean distance matrix
#'   `d_space_raw` used by the balance score.
#' @keywords internal
#' @noRd
distance_components <- function(ps, coords, E) {
  d_space_raw <- as.matrix(stats::dist(coords))

  d_space <- range01(d_space_raw)
  diag(d_space) <- 0

  D_ps <- abs(outer(ps, ps, "-"))
  D_ps <- range01(D_ps)
  diag(D_ps) <- 0

  d_exp <- as.matrix(stats::dist(E))
  d_exp <- range01(d_exp)
  diag(d_exp) <- 0

  list(D_ps = D_ps, d_space = d_space, d_exp = d_exp,
       d_space_raw = d_space_raw)
}
