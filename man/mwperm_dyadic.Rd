% Generated manually (no roxygen2/CRAN access in build environment).
\name{mwperm_dyadic}
\alias{mwperm_dyadic}
\title{Invariant permutation test for dyadic regression}
\usage{
mwperm_dyadic(y, d, x = NULL, row, col, K = NULL,
              alpha = 0.05, beta_null = 0, conf_int = TRUE,
              n_reps = 1L, seed = NULL, grid = NULL, n_cores = 1L)
}
\arguments{
\item{y}{Numeric outcome vector, one entry per observed cell.}

\item{d}{Numeric vector or matrix of the covariate(s) of interest
(\eqn{d_{ij}}). With a single covariate a confidence interval is produced;
with several, a joint confidence region (see \code{conf_int}).}

\item{x}{Optional numeric matrix or data frame of nuisance covariates
(\eqn{x_{ij}}); an intercept is always added internally. May be \code{NULL}.}

\item{row, col}{Vectors giving the row-cluster and column-cluster identity of
each observation (any type coercible by \code{\link{factor}}).}

\item{K}{Number of non-identity permutations (group order is \code{K + 1};
the smallest attainable p-value is \code{1 / (K + 1)}). Defaults to
\code{min(n_row, n_col) - 1} capped at 199, which uses the largest group
the design supports. Must satisfy \code{K + 1 <= min(n_row, n_col)}.}

\item{alpha}{Significance level used for the reject/retain decision and for
the inverted confidence interval (\code{conf_level = 1 - alpha}).}

\item{beta_null}{Null value \eqn{b} to test. Scalar (recycled if \code{d}
has several columns).}

\item{conf_int}{Logical; if \code{TRUE} a confidence set is computed by test
inversion: an interval for a single covariate, or a joint (grid-based)
confidence region for several.}

\item{n_reps}{Number of independent runs whose p-values (and CI end points)
are aggregated by the median, as recommended for randomised tests.}

\item{seed}{Optional integer; if supplied, run \code{r} uses seed
\code{seed + r - 1} for reproducibility.}

\item{grid}{Optional candidate \eqn{\beta} values for the confidence set. For a
single covariate, a numeric vector (the interval becomes the range of grid
points not rejected). For several, a list of one numeric vector per covariate
defining the region search grid.}

\item{n_cores}{Number of CPU cores for the permutation computations
(default 1 = serial). Parallelism is over the \code{n_reps} repetitions when
several are run with a \code{seed}, otherwise over the \code{K}
per-permutation factorizations; either way the result is \emph{identical} to
the serial one (every random draw is derived from explicit seeds, and the
statistics are combined by order-independent reductions). Uses forked workers
on Unix and a PSOCK cluster on Windows. Worthwhile for large problems
(thousands of cells, large \code{K}); for small ones the fork overhead
usually exceeds the gain.}
}
\value{
An object of class \code{"mwperm"} (see \code{\link{print.mwperm}}). The
provenance of its main fields: \code{estimate} and \code{se_naive} are the
\emph{OLS} point estimate(s) and the naive homoskedastic OLS standard
error(s) -- the SE is used only to centre and scale the confidence-set
search, it is not an inferential quantity; \code{conf_int} (single
coefficient) or \code{conf_region}/\code{conf_box} (several coefficients:
the retained \eqn{\beta} vectors and their marginal extent) hold the
\emph{IPT} (inverted permutation test) confidence set at level
\code{conf_level}; \code{pvalue} is the IPT permutation p-value (the median
of the per-rep p-values \code{pvalues_rep}).
}
\description{
Finite-sample valid test of \eqn{H_0: \beta = b} in the dyadic regression
model \deqn{y_{ij} = x_{ij}^\top \gamma + d_{ij}^\top \beta + \varepsilon_{ij},}
where \eqn{i} indexes the row cluster and \eqn{j} the column cluster, with a
single observation per cell. Validity holds under conditional separate
(double) exchangeability of the errors given the covariates: the
\eqn{n\times n} error matrix is row- and column-exchangeable. This is the
Invariant Permutation Test (IPT, Procedure 1) of Guo, Toulis and Wang
(2026); see \code{\link{build_perm_set}} for the permutation construction.
}
\details{
Unlike multi-way cluster-robust standard errors or the wild cluster
bootstrap, the test makes no assumption on the covariate distribution
(covariates may be irregular or heavy-tailed) and is valid for a finite
number of clusters.
}
\examples{
data(trade_dyadic)
fit <- with(trade_dyadic,
            mwperm_dyadic(y = log_trade, d = log_dist,
                          x = cbind(log_gdp_i, log_gdp_j),
                          row = importer, col = exporter,
                          seed = 1))
fit
}
\references{
Guo, F. R., Toulis, P. and Wang, Y. (2026). Permutation inference under
multi-way clustering and missing data.
}
\seealso{
\code{\link{mwperm_panel}}, \code{\link{mwperm_threeway}},
\code{\link{mwperm_layout}}, \code{\link{mwperm_missing}}.
}
