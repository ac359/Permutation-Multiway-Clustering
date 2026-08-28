#' Invariant permutation test for dyadic regression
#'
#' Finite-sample valid test of H0: beta = b in the dyadic regression
#' model
#' \preformatted{  y_ij = x_ij' gamma + d_ij' beta + eps_ij}
#'
#' where i indexes the row cluster and j the column cluster, with a single
#' observation per cell. Validity holds under conditional separate (double)
#' exchangeability of the errors given the covariates: the n-by-n error matrix
#' is row- and column-exchangeable. This is the Invariant Permutation Test (IPT,
#' Procedure 1) of Guo, Toulis and Wang (2026); see \code{\link{build_perm_set}}
#' for the permutation construction.
#'
#' Unlike multi-way cluster-robust standard errors or the wild cluster
#' bootstrap, the test makes no assumption on the covariate distribution
#' (covariates may be irregular or heavy-tailed) and is valid for a finite
#' number of clusters.
#'
#' \strong{Aggregation over repetitions, and what "exact" covers.} Theorem 1
#' establishes finite-sample validity for a \emph{single} random permutation
#' group, so with \code{n_reps = 1} the p-value is exact exactly as stated. With
#' \code{n_reps > 1} the reported p-value is the \emph{median} of the
#' per-repetition p-values, following Guo, Toulis and Wang (2026, Remark 1), and
#' the confidence set is the set of null values that the same median does not
#' reject -- one rule, applied in both places, so the test and the interval
#' cannot disagree. Each repetition's p-value is finite-sample valid on its own,
#' but no finite-sample theorem covers the median of dependent valid p-values,
#' so the median is a de-randomisation heuristic rather than a theorem.
#' Extensive simulation (every design, all of \code{n_reps} in 1--100) found the
#' median-aggregated test uniformly conservative, never anti-conservative. When
#' the guarantee must hold as stated for \code{n_reps > 1}, use \code{aggregate
#' = "median2"}: twice the median, capped at 1, is a valid p-value at level
#' \code{alpha} under arbitrary dependence across repetitions, and its
#' confidence set is correspondingly wider.
#'
#' \strong{Numerical notes.} The projection is onto the orthogonal complement of
#' the column space of the stacked design \code{[X | X_k]}, where \code{X_k} is
#' the permuted copy of \code{X}. Its dimension is therefore \code{N - rank([X |
#' X_k])}. The paper writes this as \code{N - 2p}, which is the full-rank case;
#' the stack is rank-deficient \emph{by construction} here, because a
#' permutation maps the intercept to itself, so the true dimension is strictly
#' larger than \code{N - 2p} -- substantially so in \code{\link{mwperm_panel}}
#' with \code{time_fe = TRUE}, where the period dummies are permuted among
#' themselves. Taking the rank rather than \code{2p} is the correct reading --
#' the literal \code{N - 2p} is not even well defined under rank deficiency --
#' it costs nothing, and it is what the two defining conditions on the
#' projection actually ask for. Exactly collinear nuisance columns are therefore
#' handled exactly. A \emph{nearly} collinear column of \code{x} (relative
#' tolerance about \code{1e-7} on the singular values) is dropped from the
#' projection -- the fit then behaves as if that column had been deleted, which
#' changes what is partialled out; prefer non-redundant nuisance covariates. If
#' \code{d} is annihilated entirely by some permutation's projection, the fit
#' stops rather than returning a statistic decided by rounding noise. Seeded
#' results are reproducible only under the same \code{RNGkind()} (generator and
#' sample kind) and, when cluster ids are character strings, the same collation
#' locale; see \code{\link{build_perm_set}}.
#'
#' @param y Numeric outcome vector, one entry per observed cell.
#' @param d Numeric vector or matrix of the covariate(s) of interest
#' (d_ij). With a single covariate a confidence interval is produced; with
#' several, a joint confidence region (see \code{conf_int}).
#' @param x Optional numeric matrix or data frame of nuisance covariates
#' (x_ij); an intercept is always added internally. May be \code{NULL}.
#' @param row,col Vectors giving the row-cluster and column-cluster identity of
#'   each observation (any type coercible by \code{\link{factor}}).
#' @param K Number of non-identity permutations (group order is \code{K + 1};
#' the smallest attainable p-value is \code{1 / (K + 1)}). Defaults to
#' \code{min(n_row, n_col) - 1} capped at 199, which uses the largest group the
#' design supports. Must satisfy \code{K + 1 <= min(n_row, n_col)}.
#'
#' Two things are worth knowing about that default. It maximises \emph{p-value
#' resolution}: the smallest p-value the test can return is \code{1 / (K + 1)},
#' so the largest admissible group is what makes rejection at a small
#' \code{alpha} possible at all (a 95\% confidence set needs \code{K + 1 >=
#' 20}). It also means \code{K + 1} equals the number of clusters, so Algorithm
#' 1's block decomposition is a \emph{single} block spanning every cluster, and
#' the group is one cyclic shift of the whole relabelled index set with no fixed
#' tail. The \emph{power} result, Theorem 2 of Guo, Toulis and Wang (2026), is
#' stated for a \emph{fixed} \code{K} with the number of clusters divisible by
#' \code{K + 1}, i.e. many blocks of fixed size as the design grows -- a
#' different regime. Validity is unaffected either way (Theorem 1 holds for any
#' group Algorithm 1 builds), so the default stands: it is the arrangement that
#' gives the finest grid, and the theoretical power statement simply does not
#' describe it. Set \code{K} explicitly to work in Theorem 2's regime.
#' @param alpha Significance level used for the reject/retain decision and for
#'   the inverted confidence interval (\code{conf_level = 1 - alpha}).
#' @param beta_null Null value b to test. Scalar (recycled if \code{d}
#'   has several columns).
#' @param conf_int Logical; if \code{TRUE} a confidence set is computed by test
#' inversion -- the set of null values the test does not reject, exactly as
#' Procedure 1 step 3 defines it. For a single covariate the set is computed
#' \emph{exactly} (the p-value is a step function of the null value, and the
#' package evaluates it at every jump), and its connected components are
#' returned in the \code{conf_set} field with their hull in \code{conf_int}; for
#' several covariates, a joint (grid-based) confidence region. See
#' \code{\link{confint.mwperm}}.
#' @param n_reps Number of independent runs whose p-values are aggregated by
#' \code{aggregate} (the median by default), as recommended for randomised
#' tests; the confidence set inverts the same aggregated p-value. Defaults to
#' 10: a single run's p-value depends on the random relabelling (a seed
#' lottery), and the median of 10 runs stabilises it at roughly ten times the
#' cost -- fractions of a second on typical designs. Set \code{n_reps = 1} to
#' reproduce the single-run behaviour of versions before 0.2.0. See
#' \emph{Aggregation over repetitions} in Details.
#' @param seed Optional integer; if supplied, run \code{r} uses seed
#'   \code{seed + r - 1} for reproducibility.
#' @param grid Optional candidate beta values for the confidence set. For
#' a single covariate, a numeric vector: the interval becomes the hull of the
#' grid points not rejected by the median-aggregated test (the same
#' de-randomisation as the reported p-value; an acceptance region reaching a
#' grid edge is reported as unbounded on that side). For several, a list of one
#' numeric vector per covariate defining the region search grid.
#' @param n_cores Number of CPU cores for the permutation computations
#' (default 1 = serial). Parallelism is over the \code{n_reps} repetitions when
#' several are run with a \code{seed}, otherwise over the \code{K}
#' per-permutation factorizations; either way the result is \emph{identical} to
#' the serial one (every random draw is derived from explicit seeds, and the
#' statistics are combined by order-independent reductions). Uses forked workers
#' on Unix and a PSOCK cluster on Windows. Only the seeded repetition axis
#' parallelises well (near-ideal speedup for large problems); the
#' per-permutation fallback axis is at best break-even, so with \code{n_reps =
#' 1} or \code{seed = NULL} expect little or no gain.
#'
#' @param aggregate How the \code{n_reps} per-repetition p-values are combined
#' into the reported p-value, and into the confidence set that inverts it.
#' \code{"median"} (default) is the median, as recommended in Remark 1 of Guo,
#' Toulis and Wang (2026); \code{"median2"} is \code{min(1, 2 * median)}. The
#' distinction matters for what "exact" means. Theorem 1 gives finite-sample
#' validity for a \emph{single} random permutation group, so with \code{n_reps =
#' 1} the p-value is exact as stated. The median of several dependent randomised
#' p-values is a de-randomisation heuristic -- endorsed by Remark 1 and well
#' behaved in practice, but not itself guaranteed to be a valid p-value at level
#' alpha. Twice the median is: it controls the level under arbitrary dependence
#' across repetitions (Ruschendorf 1982; Vovk and Wang 2020). Use
#' \code{"median2"} when the finite-sample guarantee must hold as stated with
#' \code{n_reps > 1}; it is conservative, never rejecting where \code{"median"}
#' would not, and its confidence set is never narrower. The default is
#' unchanged, so existing numbers stand. It costs resolution, though: since
#' \code{"median2"} reports \code{min(1, 2 * median)}, its smallest attainable
#' p-value is \code{2/(K+1)} rather than \code{1/(K+1)}, so rejecting at level
#' \code{alpha} needs \code{K + 1 >= 2/alpha} -- at alpha = 0.05 that is 40
#' levels in the smallest permuted dimension, twice the 20 \code{"median"}
#' needs. Below that the p-value is still exact but cannot reach \code{alpha},
#' and the fit reports no confidence set and says so in a note.
#' @return An object of class \code{"mwperm"} (see \code{\link{print.mwperm}}).
#' The provenance of its main fields: \code{estimate} and \code{se_naive} are
#' the \emph{OLS} point estimate(s) and the naive homoskedastic OLS standard
#' error(s) -- the SE is used only to centre and scale the confidence-set
#' search, it is not an inferential quantity; \code{conf_int} (single
#' coefficient) or \code{conf_region}/\code{conf_box} (several coefficients: the
#' retained beta vectors and their marginal extent) hold the \emph{IPT}
#' (inverted permutation test) confidence set at level \code{conf_level};
#' \code{pvalue} is the IPT permutation p-value (the median of the per-rep
#' p-values \code{pvalues_rep}).
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation
#'   inference under multi-way clustering and missing data. arXiv:2601.08610.
#'
#' @seealso \code{\link{mwperm_panel}}, \code{\link{mwperm_threeway}},
#'   \code{\link{mwperm_layout}}, \code{\link{mwperm_missing}}.
#'
#' @examples
#' data(trade_dyadic) fit <- with(trade_dyadic, mwperm_dyadic(y = log_trade, d =
#' log_dist, x = cbind(log_gdp_i, log_gdp_j), row = importer, col = exporter,
#' seed = 1)) fit
#' @export
mwperm_dyadic <- function(y, d, x = NULL, row, col, K = NULL,
                          alpha = 0.05, beta_null = 0, conf_int = TRUE,
                          n_reps = 10L, seed = NULL, grid = NULL,
                          aggregate = c("median", "median2"),
                          n_cores = 1L) {
  cl <- match.call()                   # stored on the result for printing
  aggregate <- match.arg(aggregate)
  y <- .check_y(y)
  N <- length(y)                       # number of observations
  D <- as.matrix(d)
  d_names <- .coef_names(D, deparse(substitute(d)))  # coefficient label(s)
  X <- .make_X(x, N)                   # nuisance design with intercept
  .check_lengths(N, list(row = row, col = col))

  ri <- .dense_id(row, "row")
  ci <- .dense_id(col, "col")   # dense 1-based row/col cluster ids
  n_row <- max(ri)
  n_col <- max(ci)           # number of row / col clusters
  if (anyDuplicated(cbind(ri, ci)))
    stop("Dyadic regression expects one observation per (row, col) cell. ",
         "For repeated observations use mwperm_layout() or mwperm_panel().",
         call. = FALSE)

  K <- .default_K(K, c(n_row,
                       n_col))  # group order capped by the smaller dimension
  coords <- cbind(ri, ci)              # per-observation (row, col) coordinates

  ## Per-rep permutations: draw an independent row group and column group and
  ## combine them into observation gather-vectors. Distinct sub-seeds (1, 2)
  ## keep the two dimensions' relabellings independent within a rep.
  perm_builder <- function(rep_seed) {
    Grow <- build_perm_set(n_row, K, seed = .sub_seed(rep_seed, 1L))
    Gcol <- build_perm_set(n_col, K, seed = .sub_seed(rep_seed, 2L))
    .build_obs_perms(coords, list(Grow, Gcol), design = "dyadic",
                     front_end = paste("mwperm_layout() (replication)",
                                       "or mwperm_panel() (time periods)"))
  }

  .ipt_engine(y, D, X, perm_builder, K = K, n_reps = n_reps, seed = seed,
              alpha = alpha, conf_int = conf_int, beta_null = beta_null,
              grid = grid, type = "dyadic", d_names = d_names,
              n_clusters = c(row = n_row, col = n_col), call = cl,
              n_cores = n_cores, ci_agg = aggregate)
}
