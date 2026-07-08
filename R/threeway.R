#' Invariant permutation test under three-way clustering
#'
#' Finite-sample valid test of \eqn{H_0: \beta = b} in the three-way model
#' \deqn{y_{ijl} = x_{ijl}^\top \gamma + d_{ijl}^\top \beta + \varepsilon_{ijl},}
#' \eqn{i \in [m], j \in [n], l \in [\ell]}, under full three-way
#' exchangeability of the errors (condition InvA of Guo, Toulis and Wang,
#' 2026),
#' \deqn{(\varepsilon_{ijl}) \;{\buildrel d \over =}\;
#'        (\varepsilon_{\pi(i)\sigma(j)\psi(l)}) \mid X, D,}
#' which holds, e.g., for the random-effects structure \eqn{\varepsilon_{ijl} =
#' \eta_i + \xi_j + \zeta_l + u_{ijl}} with all components i.i.d. within family.
#' A permutation is drawn independently for each of the three dimensions and
#' applied jointly.
#'
#' Use this when all three dimensions are genuinely exchangeable (e.g.
#' importer, exporter and product category). If one dimension is time with
#' autocorrelation, use \code{\link{mwperm_panel}} instead. The data must form
#' a complete balanced array.
#'
#' @inheritParams mwperm_dyadic
#' @param id1,id2,id3 Cluster identifiers for the three dimensions.
#' @param K Number of non-identity permutations; defaults to
#'   \code{min(m, n, ell) - 1} capped at 199.
#'
#' @return An object of class \code{"mwperm"}: \code{estimate}/\code{se_naive}
#'   are the OLS estimate and naive SE, \code{conf_int} (or
#'   \code{conf_region}/\code{conf_box} for several coefficients) the IPT
#'   inverted-test confidence set, and \code{pvalue} the IPT permutation
#'   p-value; see \code{\link{mwperm_dyadic}} for the field provenance in full.
#' @references Guo, F. R., Toulis, P. and Wang, Y. (2026), Section 6.1.
#' @seealso \code{\link{mwperm_dyadic}}, \code{\link{mwperm_panel}}.
#' @export
mwperm_threeway <- function(y, d, x = NULL, id1, id2, id3, K = NULL,
                            alpha = 0.05, beta_null = 0, conf_int = TRUE,
                            n_reps = 1L, seed = NULL, grid = NULL,
                            n_cores = 1L) {
  cl <- match.call()
  y <- as.numeric(y); N <- length(y)
  D <- as.matrix(d); d_names <- .coef_names(D, deparse(substitute(d)))
  X <- .make_X(x, N)
  .check_lengths(N, list(id1 = id1, id2 = id2, id3 = id3))

  a1 <- .dense_id(id1, "id1"); a2 <- .dense_id(id2, "id2")          # dense ids per dimension
  a3 <- .dense_id(id3, "id3")
  m <- max(a1); n <- max(a2); ell <- max(a3)                        # cluster counts per dimension
  coords <- cbind(a1, a2, a3)          # per-observation (id1, id2, id3) coordinates
  .require_complete_array(coords, c(id1 = m, id2 = n, id3 = ell), N,
                          what = "Three-way design")

  K <- .default_K(K, c(m, n, ell))     # group order capped by the smallest dimension

  ## An independent permutation is drawn for each of the three dimensions (distinct
  ## sub-seeds 1, 2, 3) and applied jointly: full three-way exchangeability (InvA).
  perm_builder <- function(rep_seed) {
    G1 <- build_perm_set(m,   K, seed = .sub_seed(rep_seed, 1L))
    G2 <- build_perm_set(n,   K, seed = .sub_seed(rep_seed, 2L))
    G3 <- build_perm_set(ell, K, seed = .sub_seed(rep_seed, 3L))
    .build_obs_perms(coords, list(G1, G2, G3))
  }

  .ipt_engine(y, D, X, perm_builder, K = K, n_reps = n_reps, seed = seed,
              alpha = alpha, conf_int = conf_int, beta_null = beta_null,
              grid = grid, type = "threeway", d_names = d_names,
              n_clusters = c(id1 = m, id2 = n, id3 = ell), call = cl,
              n_cores = n_cores)
}
