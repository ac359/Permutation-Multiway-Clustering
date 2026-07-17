#' Invariant permutation test for panel (longitudinal) dyadic regression
#'
#' Finite-sample valid test of \eqn{H_0: \beta = b} in the panel model
#' \deqn{y_{ijt} = x_{ijt}^\top \gamma + d_{ijt}^\top \beta + \varepsilon_{ijt},}
#' where \eqn{i, j} index two cross-sectional clustering dimensions (e.g.
#' importer and exporter countries) and \eqn{t} indexes time. Full three-way
#' exchangeability is implausible because errors are typically autocorrelated
#' over time. Instead this test only assumes exchangeability across the first
#' two dimensions \emph{within} each time period (condition InvB of Guo, Toulis
#' and Wang, 2026),
#' \deqn{(\varepsilon_{ijt})_{i,j} \;{\buildrel d \over =}\;
#'        (\varepsilon_{\pi(i)\sigma(j)t})_{i,j} \mid X, D,}
#' which holds, for instance, under \eqn{\varepsilon_{ijt} = \eta_i + \xi_j +
#' \zeta_t + u_{ijt}} with an \emph{arbitrary} time trend \eqn{\zeta_t}. The
#' same row/column permutation is applied in every period, so any unknown time
#' effect is held fixed and partialled out. This is, to the authors' knowledge,
#' the first finite-sample valid test for \eqn{\beta = 0} under exchangeable
#' errors in such panel models.
#'
#' The data must form a complete balanced array: every \eqn{(i, j, t)} cell
#' present exactly once.
#'
#' A feasibility note for long panels: the engine requires \code{N > 2p}
#' (with \eqn{p} the number of nuisance columns including the intercept and,
#' with \code{time_fe = TRUE}, the time dummies). This matches the premise of
#' the validity theorem and is slightly conservative here -- because time is
#' held fixed, the time dummies are duplicated between the design and its
#' permuted copy, so the stacked projection's true rank is below \eqn{2p}.
#'
#' @inheritParams mwperm_dyadic
#' @param d Numeric vector or matrix of the covariate(s) of interest
#'   \eqn{d_{ijt}} (may be time-varying). With a single covariate a confidence
#'   interval is produced; with several, a joint confidence region.
#' @param x Optional nuisance covariates \eqn{x_{ijt}}; intercept added
#'   internally. May be \code{NULL}.
#' @param row,col Row- and column-cluster identifiers (length = number of
#'   observations).
#' @param time Time-period identifiers (length = number of observations).
#' @param time_fe Logical; if \code{TRUE} (the default) time fixed effects are
#'   added to the nuisance design. Under condition InvB this is valid (the time
#'   dummies are invariant to the within-period permutation) and it removes the
#'   time trend \eqn{\zeta_t} from the residuals, which de-biases the reported
#'   point estimate and sharpens the test when treatment timing is correlated
#'   with the period.
#' @param K Number of non-identity permutations; defaults to
#'   \code{min(n_row, n_col) - 1} capped at 199.
#'
#' @return An object of class \code{"mwperm"}: \code{estimate}/\code{se_naive}
#'   are the OLS estimate and naive SE, \code{conf_int} (or
#'   \code{conf_region}/\code{conf_box} for several coefficients) the IPT
#'   inverted-test confidence set, and \code{pvalue} the IPT permutation
#'   p-value; see \code{\link{mwperm_dyadic}} for the field provenance in full.
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation
#'   inference under multi-way clustering and missing data, Section 6.2.
#'   arXiv:2601.08610.
#'
#' @seealso \code{\link{mwperm_dyadic}}, \code{\link{mwperm_threeway}}.
#'
#' @examples
#' data(trade_panel)
#' fit <- with(trade_panel,
#'             mwperm_panel(y = log_trade, d = fta, x = cbind(log_gdp_i, log_gdp_j),
#'                          row = importer, col = exporter, time = year, seed = 1))
#' fit
#' @export
mwperm_panel <- function(y, d, x = NULL, row, col, time, K = NULL,
                         alpha = 0.05, beta_null = 0, conf_int = TRUE,
                         n_reps = 10L, seed = NULL, grid = NULL, time_fe = TRUE,
                         n_cores = 1L) {
  cl <- match.call()
  y <- .check_y(y); N <- length(y)
  D <- as.matrix(d); d_names <- .coef_names(D, deparse(substitute(d)))
  .check_lengths(N, list(row = row, col = col, time = time))

  ri <- .dense_id(row, "row"); ci <- .dense_id(col, "col")           # dense ids per dimension
  ti <- .dense_id(time, "time")
  n_row <- max(ri); n_col <- max(ci); n_t <- max(ti)                 # cluster / period counts

  ## Nuisance design: user covariates (+ intercept) and optional time dummies.
  ## The time dummies are invariant to the within-period permutation, so adding
  ## them is valid under condition InvB and removes the time trend zeta_t.
  X <- .make_X(x, N)
  if (isTRUE(time_fe) && n_t > 1L) {
    TD <- stats::model.matrix(~ factor(ti))[, -1L, drop = FALSE]   # period dummies (drop reference)
    colnames(TD) <- paste0("time", sort(unique(ti))[-1L])
    X <- cbind(X, TD)
  }

  coords <- cbind(ri, ci, ti)
  .require_complete_array(coords, c(row = n_row, col = n_col, time = n_t), N,
                          what = "Panel")

  K <- .default_K(K, c(n_row, n_col))

  ## The SAME row/column permutation is applied in every period (time passed as
  ## NULL = held fixed), so any unknown time effect is preserved and partialled out.
  perm_builder <- function(rep_seed) {
    Grow <- build_perm_set(n_row, K, seed = .sub_seed(rep_seed, 1L))
    Gcol <- build_perm_set(n_col, K, seed = .sub_seed(rep_seed, 2L))
    .build_obs_perms(coords, list(Grow, Gcol, NULL))   # time held fixed
  }

  .ipt_engine(y, D, X, perm_builder, K = K, n_reps = n_reps, seed = seed,
              alpha = alpha, conf_int = conf_int, beta_null = beta_null,
              grid = grid, type = "panel", d_names = d_names,
              n_clusters = c(row = n_row, col = n_col, time = n_t), call = cl,
              n_cores = n_cores)
}
