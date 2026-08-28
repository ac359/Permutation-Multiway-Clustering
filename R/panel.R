#' Invariant permutation test for panel (longitudinal) dyadic regression
#'
#' Finite-sample valid test of H0: beta = b in the panel model
#' \preformatted{  y_ijt = x_ijt' gamma + d_ijt' beta + eps_ijt}
#'
#' where i, j index two cross-sectional clustering dimensions (e.g.
#' importer and exporter countries) and t indexes time. Full three-way
#' exchangeability is implausible because errors are typically autocorrelated
#' over time. Instead this test only assumes exchangeability across the first
#' two dimensions \emph{within} each time period (condition InvB of Guo, Toulis
#' and Wang, 2026),
#' \preformatted{  (eps_ijt) has the same distribution as
#'     (eps_pi(i),sigma(j),t), for every row permutation pi and column
#'     permutation sigma, conditional on X and D, with t held fixed}
#'
#' which holds, for instance, under eps_ijt = eta_i + xi_j + zeta_t + u_ijt with
#' an \emph{arbitrary} time trend zeta_t. The same row/column permutation is
#' applied in every period, so any unknown time effect is held fixed and
#' partialled out. This is, to the authors' knowledge, the first finite-sample
#' valid test for beta = 0 under exchangeable errors in such panel models.
#'
#' The data must form a complete balanced array: every (i, j, t) cell present
#' exactly once.
#'
#' A feasibility note for long panels: the engine requires \code{N > 2p} (with p
#' the number of nuisance columns including the intercept and, with
#' \code{time_fe = TRUE}, the time dummies). This matches the premise of the
#' validity theorem and is conservative here. The projection is onto the
#' orthogonal complement of the column space of \code{[X | X_k]}, where
#' \code{X_k} is the permuted copy of \code{X}, so its dimension is \code{N -
#' rank([X | X_k])}, not the \code{N - 2p} of the paper's full-rank statement:
#' the intercept is permutation-invariant in every design, and here the time
#' dummies are too (time is held fixed), so the stack is rank-deficient by
#' construction and the projection retains \emph{more} dimensions than \code{N -
#' 2p}. That is the correct reading -- \code{N - 2p} is not well defined when
#' the stack is rank deficient -- and it is why \code{time_fe = TRUE} costs far
#' less than its column count suggests.
#'
#' @inheritParams mwperm_dyadic
#' @param d Numeric vector or matrix of the covariate(s) of interest
#' d_ijt (may be time-varying). With a single covariate a confidence interval is
#' produced; with several, a joint confidence region.
#' @param x Optional nuisance covariates x_ijt; intercept added
#'   internally. May be \code{NULL}.
#' @param row,col Row- and column-cluster identifiers (length = number of
#'   observations).
#' @param time Time-period identifiers (length = number of observations).
#' @param time_fe Logical; if \code{TRUE} (the default) time fixed effects are
#' added to the nuisance design. Under condition InvB this is valid (the time
#' dummies are invariant to the within-period permutation) and it removes the
#' time trend zeta_t from the residuals, which de-biases the reported point
#' estimate and sharpens the test when treatment timing is correlated with the
#' period.
#' @param K Number of non-identity permutations; defaults to
#'   \code{min(n_row, n_col) - 1} capped at 199.
#'
#' @param aggregate How the \code{n_reps} per-repetition p-values are
#'   combined into the reported p-value, and into the confidence set that
#'   inverts it. \code{"median"} (the default) is the median, as recommended in
#'   Remark 1 of Guo, Toulis and Wang (2026); \code{"median2"} is
#'   \code{min(1, 2 * median)}.
#'
#'   The choice decides what "exact" covers. Theorem 1 gives finite-sample
#'   validity for a single random permutation group, so at \code{n_reps = 1}
#'   the p-value is exact as stated. The median of several dependent randomised
#'   p-values is a de-randomisation heuristic: endorsed by Remark 1 and well
#'   behaved in practice, but not itself guaranteed valid at level
#'   \code{alpha}. Twice the median is guaranteed, under arbitrary dependence
#'   across repetitions (Ruschendorf 1982; Vovk and Wang 2020).
#'
#'   So use \code{"median2"} when the guarantee must hold as stated with
#'   \code{n_reps > 1}. It is conservative: it never rejects where
#'   \code{"median"} would not, and its confidence set is never narrower. The
#'   default is unchanged, so existing numbers stand.
#'
#'   The cost is resolution. \code{"median2"} reports
#'   \code{min(1, 2 * median)}, so its smallest attainable p-value is
#'   \code{2/(K+1)}, not \code{1/(K+1)}, and rejecting at level \code{alpha}
#'   needs \code{K + 1 >= 2/alpha} -- at \code{alpha = 0.05} that is 40 levels
#'   in the smallest permuted dimension, twice what \code{"median"} needs.
#'   Below that the p-value is still exact but cannot reach \code{alpha}, and
#'   the fit says so in a note.
#' @return An object of class \code{"mwperm"}: \code{estimate}/\code{se_naive}
#' are the OLS estimate and naive SE, \code{conf_int} (or
#' \code{conf_region}/\code{conf_box} for several coefficients) the IPT
#' inverted-test confidence set, and \code{pvalue} the IPT permutation p-value;
#' see \code{\link{mwperm_dyadic}} for the field provenance in full.
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation
#' inference under multi-way clustering and missing data, Section 6.2.
#' arXiv:2601.08610.
#'
#' @seealso \code{\link{mwperm_dyadic}}, \code{\link{mwperm_threeway}}.
#'
#' @examples
#' data(trade_panel) fit <- with(trade_panel, mwperm_panel(y = log_trade, d =
#' fta, x = cbind(log_gdp_i, log_gdp_j), row = importer, col = exporter, time =
#' year, seed = 1)) fit
#' @export
mwperm_panel <- function(y, d, x = NULL, row, col, time, K = NULL,
                         alpha = 0.05, beta_null = 0, conf_int = TRUE,
                         n_reps = 10L, seed = NULL, grid = NULL, time_fe = TRUE,
                         aggregate = c("median", "median2"),
                         n_cores = 1L) {
  cl <- match.call()
  aggregate <- match.arg(aggregate)
  y <- .check_y(y)
  N <- length(y)
  D <- as.matrix(d)
  d_names <- .coef_names(D, deparse(substitute(d)))
  .check_lengths(N, list(row = row, col = col, time = time))

  ri <- .dense_id(row, "row")
  ci <- .dense_id(col, "col")           # dense ids per dimension
  ti <- .dense_id(time, "time")
  n_row <- max(ri)
  n_col <- max(ci)
  n_t <- max(ti)                 # cluster / period counts

  ## Nuisance design: user covariates (+ intercept) and optional time dummies.
  ## The time dummies are invariant to the within-period permutation, so adding
  ## them is valid under condition InvB and removes the time trend zeta_t.
  X <- .make_X(x, N)
  if (isTRUE(time_fe) && n_t > 1L) {
    ## period dummies (drop the reference level)
    TD <- stats::model.matrix(~ factor(ti))[, -1L, drop = FALSE]
    colnames(TD) <- paste0("time", sort(unique(ti))[-1L])
    X <- cbind(X, TD)
  }

  coords <- cbind(ri, ci, ti)
  .require_complete_array(coords, c(row = n_row, col = n_col, time = n_t), N,
                          what = "Panel")

  K <- .default_K(K, c(n_row, n_col))

  ## The SAME row/column permutation is applied in every period (time passed as
  ## NULL = held fixed), so any unknown time effect is preserved and partialled
  ## out.
  perm_builder <- function(rep_seed) {
    Grow <- build_perm_set(n_row, K, seed = .sub_seed(rep_seed, 1L))
    Gcol <- build_perm_set(n_col, K, seed = .sub_seed(rep_seed, 2L))
    ## time held fixed (the NULL group) -- this is condition InvB
    .build_obs_perms(coords, list(Grow, Gcol, NULL), design = "panel",
                     front_end = "mwperm_layout() (within-cell replication)")
  }

  .ipt_engine(y, D, X, perm_builder, K = K, n_reps = n_reps, seed = seed,
              alpha = alpha, conf_int = conf_int, beta_null = beta_null,
              grid = grid, type = "panel", d_names = d_names,
              n_clusters = c(row = n_row, col = n_col, time = n_t), call = cl,
              n_cores = n_cores, ci_agg = aggregate)
}
