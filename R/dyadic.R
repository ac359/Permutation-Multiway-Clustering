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
#' @param seed Optional integer; if supplied, run \code{r} uses seed
#'   \code{seed + r - 1} for reproducibility.
#' @param grid Optional candidate \eqn{\beta} values for the confidence set. For
#'   a single covariate, a numeric vector (the interval becomes the range of grid
#'   points not rejected). For several, a list of one numeric vector per
#'   covariate defining the region search grid.
#'
#' @return An object of class \code{"mwperm"} (see \code{\link{print.mwperm}}).
#'
#' @references Guo, F. R., Toulis, P. and Wang, Y. (2026). Permutation
#'   inference under multi-way clustering and missing data.
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
                          n_reps = 1L, seed = NULL, grid = NULL) {
  cl <- match.call()
  y <- as.numeric(y)
  N <- length(y)
  D <- as.matrix(d); d_names <- .coef_names(D, deparse(substitute(d)))
  X <- .make_X(x, N)
  .check_lengths(N, list(row = row, col = col))

  ri <- .dense_id(row); ci <- .dense_id(col)
  n_row <- max(ri); n_col <- max(ci)
  if (anyDuplicated(cbind(ri, ci)))
    stop("Dyadic regression expects one observation per (row, col) cell. ",
         "For repeated observations use mwperm_layout() or mwperm_panel().",
         call. = FALSE)

  K <- .default_K(K, c(n_row, n_col))
  coords <- cbind(ri, ci)

  perm_builder <- function(rep_seed) {
    Grow <- build_perm_set(n_row, K, seed = .sub_seed(rep_seed, 1L))
    Gcol <- build_perm_set(n_col, K, seed = .sub_seed(rep_seed, 2L))
    .build_obs_perms(coords, list(Grow, Gcol))
  }

  .ipt_engine(y, D, X, perm_builder, K = K, n_reps = n_reps, seed = seed,
              alpha = alpha, conf_int = conf_int, beta_null = beta_null,
              grid = grid, conf_level = 1 - alpha, type = "dyadic",
              d_names = d_names,
              n_clusters = c(row = n_row, col = n_col), call = cl)
}

## ---- small shared helpers -------------------------------------------------

.default_K <- function(K, dim_sizes, cap = 199L) {
  gmax <- min(dim_sizes) - 1L
  if (is.null(K)) {
    K <- min(gmax, cap)
  } else {
    K <- as.integer(K)
    if (K + 1L > min(dim_sizes))
      stop(sprintf("K + 1 = %d exceeds the smallest permuted dimension (%d).",
                   K + 1L, min(dim_sizes)), call. = FALSE)
  }
  if (K < 1L) stop("Not enough clusters to permute (need >= 2 in each dimension).",
                   call. = FALSE)
  K
}

.sub_seed <- function(rep_seed, j) if (is.null(rep_seed)) NULL else rep_seed * 1000L + j

.coef_names <- function(D, fallback) {
  nm <- colnames(D)
  if (!is.null(nm)) return(nm)
  if (ncol(D) == 1L) fallback else paste0(fallback, seq_len(ncol(D)))
}
