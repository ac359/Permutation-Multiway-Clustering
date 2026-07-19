#' Invariant permutation test for dyadic regression
#'
#' Finite-sample valid test of \eqn{H_0: \beta = b} in the dyadic regression
#' model
#' \deqn{y_{ij} = x_{ij}^\top \gamma + d_{ij}^\top \beta + \varepsilon_{ij},}
#' where \eqn{i} indexes the row cluster and \eqn{j} the column cluster, with a
#' single observation per cell. Validity holds under conditional separate
#' (double) exchangeability of the errors given the covariates: the
#' \eqn{n\times n} error matrix is row- and column-exchangeable. This is the
#' Invariant Permutation Test (IPT, Procedure 1) of Guo, Toulis and Wang
#' (2026); see \code{\link{build_perm_set}} for the permutation construction.
#'
#' Unlike multi-way cluster-robust standard errors or the wild cluster
#' bootstrap, the test makes no assumption on the covariate distribution
#' (covariates may be irregular or heavy-tailed) and is valid for a finite
#' number of clusters.
#'
#' \strong{Aggregation over repetitions.} With \code{n_reps > 1} the reported
#' p-value is the \emph{median} of the per-repetition p-values and the
#' confidence limits are medians of the per-repetition endpoints, following
#' Guo, Toulis and Wang (2026, Remark 1). Each repetition's p-value is
#' finite-sample valid on its own; no finite-sample theorem covers the median
#' of dependent valid p-values in general, so the aggregated p-value should
#' be read as a stabilised summary. Extensive simulation (every design, all
#' of \code{n_reps} in 1--100) found the median-aggregated test uniformly
#' conservative, never anti-conservative; a provably valid fallback is to
#' reject only when \emph{twice} the median p-value is at most \code{alpha}
#' (see the aggregation references cited in Remark 1 of the paper).
#'
#' \strong{Numerical notes.} Exactly collinear nuisance columns (such as the
#' intercept duplicated in the stacked projection) are handled exactly by the
#' rank-revealing QR. A \emph{nearly} collinear column of \code{x} (relative
#' QR tolerance about \code{1e-7}) is silently dropped from the projection --
#' the fit then behaves exactly as if that column had been deleted, which
#' changes what is partialled out; prefer non-redundant nuisance covariates.
#' Seeded results are reproducible only under the same \code{RNGkind()}
#' (generator and sample kind) and, when cluster ids are character strings,
#' the same collation locale; see \code{\link{build_perm_set}}.
#'
#' @param y Numeric outcome vector, one entry per observed cell.
#' @param d Numeric vector or matrix of the covariate(s) of interest
#'   (\eqn{d_{ij}}). With a single covariate a confidence interval is produced;
#'   with several, a joint confidence region (see \code{conf_int}).
#' @param x Optional numeric matrix or data frame of nuisance covariates
#'   (\eqn{x_{ij}}); an intercept is always added internally. May be
#'   \code{NULL}.
#' @param row,col Vectors giving the row-cluster and column-cluster identity of
#'   each observation (any type coercible by \code{\link{factor}}).
#' @param K Number of non-identity permutations (group order is \code{K + 1};
#'   the smallest attainable p-value is \code{1 / (K + 1)}). Defaults to
#'   \code{min(n_row, n_col) - 1} capped at 199, which uses the largest group
#'   the design supports. Must satisfy \code{K + 1 <= min(n_row, n_col)}.
#' @param alpha Significance level used for the reject/retain decision and for
#'   the inverted confidence interval (\code{conf_level = 1 - alpha}).
#' @param beta_null Null value \eqn{b} to test. Scalar (recycled if \code{d}
#'   has several columns).
#' @param conf_int Logical; if \code{TRUE} a confidence set is computed by test
#'   inversion: an interval for a single covariate, or a joint (grid-based)
#'   confidence region for several.
#' @param n_reps Number of independent runs whose p-values (and CI end points)
#'   are aggregated by the median, as recommended for randomised tests.
#'   Defaults to 10: a single run's p-value depends on the random relabelling
#'   (a seed lottery), and the median of 10 runs stabilises it at roughly ten
#'   times the cost -- fractions of a second on typical designs. Set
#'   \code{n_reps = 1} to reproduce the single-run behaviour of versions
#'   before 0.2.0. See \emph{Aggregation over repetitions} in Details.
#' @param seed Optional integer; if supplied, run \code{r} uses seed
#'   \code{seed + r - 1} for reproducibility.
#' @param grid Optional candidate \eqn{\beta} values for the confidence set. For
#'   a single covariate, a numeric vector: the interval becomes the hull of the
#'   grid points not rejected by the median-aggregated test (the same
#'   de-randomisation as the reported p-value; an acceptance region reaching a
#'   grid edge is reported as unbounded on that side). For several, a list of
#'   one numeric vector per covariate defining the region search grid.
#' @param n_cores Number of CPU cores for the permutation computations
#'   (default 1 = serial). Parallelism is over the \code{n_reps} repetitions
#'   when several are run with a \code{seed}, otherwise over the \code{K}
#'   per-permutation factorizations; either way the result is \emph{identical}
#'   to the serial one (every random draw is derived from explicit seeds, and
#'   the statistics are combined by order-independent reductions). Uses forked
#'   workers on Unix and a PSOCK cluster on Windows. Only the seeded
#'   repetition axis parallelises well (near-ideal speedup for large
#'   problems); the per-permutation fallback axis is at best break-even, so
#'   with \code{n_reps = 1} or \code{seed = NULL} expect little or no gain.
#'
#' @return An object of class \code{"mwperm"} (see \code{\link{print.mwperm}}).
#'   The provenance of its main fields: \code{estimate} and \code{se_naive}
#'   are the \emph{OLS} point estimate(s) and the naive homoskedastic OLS
#'   standard error(s) -- the SE is used only to centre and scale the
#'   confidence-set search, it is not an inferential quantity;
#'   \code{conf_int} (single coefficient) or \code{conf_region}/\code{conf_box}
#'   (several coefficients: the retained \eqn{\beta} vectors and their
#'   marginal extent) hold the \emph{IPT} (inverted permutation test)
#'   confidence set at level \code{conf_level}; \code{pvalue} is the IPT
#'   permutation p-value (the median of the per-rep p-values
#'   \code{pvalues_rep}).
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation
#'   inference under multi-way clustering and missing data. arXiv:2601.08610.
#'
#' @seealso \code{\link{mwperm_panel}}, \code{\link{mwperm_threeway}},
#'   \code{\link{mwperm_layout}}, \code{\link{mwperm_missing}}.
#'
#' @examples
#' data(trade_dyadic)
#' fit <- with(trade_dyadic,
#'             mwperm_dyadic(y = log_trade, d = log_dist,
#'                           x = cbind(log_gdp_i, log_gdp_j),
#'                           row = importer, col = exporter,
#'                           seed = 1))
#' fit
#' @export
mwperm_dyadic <- function(y, d, x = NULL, row, col, K = NULL,
                          alpha = 0.05, beta_null = 0, conf_int = TRUE,
                          n_reps = 10L, seed = NULL, grid = NULL,
                          n_cores = 1L) {
  cl <- match.call()                   # stored on the result for printing
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
    .build_obs_perms(coords, list(Grow, Gcol))
  }

  .ipt_engine(y, D, X, perm_builder, K = K, n_reps = n_reps, seed = seed,
              alpha = alpha, conf_int = conf_int, beta_null = beta_null,
              grid = grid, type = "dyadic", d_names = d_names,
              n_clusters = c(row = n_row, col = n_col), call = cl,
              n_cores = n_cores)
}
