## ============================================================================
## R/threeway.R -- the three-way design: Section 6.1, condition InvA
##
## Purpose. mwperm_threeway() validates a complete m x n x ell array and
##   hands the engine a builder of three independent Algorithm 1 groups
##   applied jointly.
## Paper. Guo, Toulis & Wang (2026), Section 6.1: model (12), condition
##   InvA (full three-way exchangeability, e.g. three-way random effects);
##   step (i) Algorithm 1 three times, step (ii) Procedure 1. Not for a time
##   index with autocorrelation -- that is Section 6.2 (panel.R), and
##   mwperm() never picks this design unless asked.
## Pipeline. mwperm() -> dispatch -> [design worker] -> permutation
##   construction -> projection engine -> median aggregation -> test
##   inversion -> S3 methods.
## ============================================================================

#' Invariant permutation test under three-way clustering
#'
#' Finite-sample valid test of H0: beta = b in the three-way model
#' \preformatted{  y_ijl = x_ijl' gamma + d_ijl' beta + eps_ijl}
#'
#' with i in `[m]`, j in `[n]` and l in `[ell]`, under full three-way
#' exchangeability of the errors (condition InvA of Guo, Toulis and Wang,
#' 2026),
#' \preformatted{  (eps_ijl) has the same distribution as
#'     (eps_pi(i),sigma(j),psi(l)), for every permutation pi of the rows,
#'     sigma of the columns and psi of the third dimension, given X and D}
#'
#' which holds, e.g., for the random-effects structure eps_ijl = eta_i + xi_j
#' + zeta_l + u_ijl with all components i.i.d. within family. A permutation is
#' drawn independently for each of the three dimensions and applied jointly.
#'
#' Use this when all three dimensions are genuinely exchangeable (e.g.
#' importer, exporter and product category). If one dimension is time with
#' autocorrelation, use [mwperm_panel()] instead. The data must form a
#' complete balanced array.
#'
#' @details Implements Section 6.1 of Guo, Toulis and Wang (2026): Algorithm 1
#'   applied independently to each of the three dimensions, the elements
#'   applied jointly (condition InvA), then Procedure 1.
#' @inheritParams mwperm_dyadic
#' @param id1,id2,id3 Cluster identifiers for the three dimensions.
#' @param K Number of non-identity permutations; defaults to `min(m, n, ell) -
#'   1` capped at 199.
#'
#' @param aggregate How the `n_reps` per-repetition p-values are combined into
#'   the reported p-value, and into the confidence set that inverts it.
#'   `"median"` (the default) is the median, as recommended in Remark 1 of
#'   Guo, Toulis and Wang (2026); `"median2"` is `min(1, 2 * median)`.
#'
#' The choice decides what "exact" covers. Theorem 1 gives finite-sample
#' validity for a single random permutation group, so at `n_reps = 1` the
#' p-value is exact as stated. The median of several dependent randomised
#' p-values is a de-randomisation heuristic: endorsed by Remark 1 and well
#' behaved in practice, but not itself guaranteed valid at level `alpha`.
#' Twice the median is guaranteed, under arbitrary dependence across
#' repetitions (Ruschendorf 1982; Vovk and Wang 2020).
#'
#' So use `"median2"` when the guarantee must hold as stated with `n_reps >
#' 1`. It is conservative: it never rejects where `"median"` would not, and
#' its confidence set is never narrower. The default is unchanged, so existing
#' numbers stand.
#'
#' The cost is resolution. `"median2"` reports `min(1, 2 * median)`, so its
#' smallest attainable p-value is `2/(K+1)`, not `1/(K+1)`, and rejecting at
#' level `alpha` needs `K + 1 >= 2/alpha` -- at `alpha = 0.05` that is 40
#' levels in the smallest permuted dimension, twice what `"median"` needs.
#' Below that the p-value is still exact but cannot reach `alpha`, and the fit
#' says so in a note.
#' @return An object of class `"mwperm"`: `estimate`/`se_naive` are the OLS
#'   estimate and naive SE, `conf_int` (or `conf_region`/`conf_box` for
#'   several coefficients) the IPT inverted-test confidence set, and `pvalue`
#'   the IPT permutation p-value; see [mwperm_dyadic()] for the field
#'   provenance in full.
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data, Section 6.1.
#'   arXiv:2601.08610.
#' @seealso [mwperm_dyadic()], [mwperm_panel()].
#' @examples
#' ## small balanced 3-way array with a random-effects error
#' set.seed(1)
#' m <- 8
#' g <- expand.grid(i = seq_len(m), j = seq_len(m), l = seq_len(m))
#' a <- rnorm(m); b <- rnorm(m); cc <- rnorm(m)
#' g$d <- rnorm(nrow(g))
#' g$y <- 0.6 * g$d + a[g$i] + b[g$j] + cc[g$l] + rnorm(nrow(g))
#' fit <- with(g, mwperm_threeway(y = y, d = d, id1 = i, id2 = j, id3 = l,
#'                                conf_int = FALSE, seed = 1))
#' fit
#' @export
mwperm_threeway <- function(y, d, x = NULL, id1, id2, id3, K = NULL,
                            alpha = 0.05, beta_null = 0, conf_int = TRUE,
                            n_reps = 10L, seed = NULL, grid = NULL,
                            aggregate = c("median", "median2"),
                            n_cores = 1L) {
  cl <- match.call()
  aggregate <- match.arg(aggregate)
  y <- .check_y(y)
  N <- length(y)
  D <- as.matrix(d)
  d_names <- .coef_names(D, deparse(substitute(d)))
  X <- .make_X(x, N)
  .check_lengths(N, list(id1 = id1, id2 = id2, id3 = id3))

  a1 <- .dense_id(id1, "id1")
  a2 <- .dense_id(id2, "id2")          # dense ids per dimension
  a3 <- .dense_id(id3, "id3")
  m <- max(a1)
  n <- max(a2)
  ell <- max(a3)                        # cluster counts per dimension
  coords <- cbind(a1, a2,
                  a3)          # per-observation (id1, id2, id3) coordinates
  .require_complete_array(coords, c(id1 = m, id2 = n, id3 = ell), N,
                          what = "Three-way design")

  K <- .default_K(K, c(m, n,
                       ell))     # group order capped by the smallest dimension

  ## An independent permutation is drawn for each of the three dimensions
  ## (distinct
  ## sub-seeds 1, 2, 3) and applied jointly: full three-way exchangeability
  ## (InvA).
  perm_builder <- function(rep_seed) {
    ## Section 6.1, step (i): Algorithm 1 applied three times, one group per
    ## dimension; element k acts as (psi^1_k, psi^2_k, psi^3_k) jointly.
    G1 <- build_perm_set(m,   K, seed = .sub_seed(rep_seed, 1L))
    G2 <- build_perm_set(n,   K, seed = .sub_seed(rep_seed, 2L))
    G3 <- build_perm_set(ell, K, seed = .sub_seed(rep_seed, 3L))
    .build_obs_perms(coords, list(G1, G2, G3), design = "threeway",
                     front_end = "mwperm_layout() (within-cell replication)")
  }

  .ipt_engine(y, D, X, perm_builder, K = K, n_reps = n_reps, seed = seed,
              alpha = alpha, conf_int = conf_int, beta_null = beta_null,
              grid = grid, type = "threeway", d_names = d_names,
              n_clusters = c(id1 = m, id2 = n, id3 = ell), call = cl,
              n_cores = n_cores, ci_agg = aggregate)
}
