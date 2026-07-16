% Generated manually (no roxygen2/CRAN access in build environment).
\name{mwperm_missing}
\alias{mwperm_missing}
\title{Invariant permutation test for dyadic regression with missing cells}
\usage{
mwperm_missing(y, d, x = NULL, row, col, K = NULL,
               alpha = 0.05, beta_null = 0, conf_int = TRUE,
               n_reps = 1L, seed = NULL, grid = NULL, min_block = 3L,
               block_method = c("greedy", "exact"), n_cores = 1L)
}
\arguments{
\item{y}{Numeric outcome vector, one entry per observed cell.}

\item{d}{Numeric vector or matrix of the covariate(s) of interest. With a
single covariate a confidence interval is produced; with several, a joint
confidence region.}

\item{x}{Optional numeric matrix or data frame of nuisance covariates; an
intercept is always added internally. May be \code{NULL}.}

\item{row, col}{Cluster identities of each observed cell. Cells that are
absent from the data are treated as missing.}

\item{K}{Number of non-identity permutations. Defaults to
\code{min_q min(|I_q|, |J_q|) - 1} over the selected blocks, capped at 199.
Must not exceed that quantity.}

\item{alpha}{Significance level for the decision and the inverted interval.}

\item{beta_null}{Null value \eqn{b} to test.}

\item{conf_int}{Logical; compute a test-inversion confidence set: an interval
for a single covariate, or a joint confidence region for several.}

\item{n_reps}{Number of independent runs aggregated by the median.}

\item{seed}{Optional integer; run \code{r} uses seed \code{seed + r - 1}.}

\item{grid}{Optional numeric vector of candidate \eqn{\beta} values.}
\item{n_cores}{Number of CPU cores for the permutation computations
(default 1 = serial). Parallelism is over the \code{n_reps} repetitions when
several are run with a \code{seed}, otherwise over the \code{K}
per-permutation factorizations; either way the result is \emph{identical} to
the serial one. See \code{\link{mwperm_dyadic}} for details.}

\item{min_block}{Integer; only bicliques whose smaller side is at least
\code{min_block} are used. Larger values give finer p-value resolution
(larger \code{K}) but discard more data. Defaults to 3.}

\item{block_method}{Passed to \code{\link{find_bicliques}}: \code{"greedy"}
(default) for the fast heuristic, or \code{"exact"} for a branch-and-bound
search of maximum-area blocks (with a node-budget fallback to greedy). See
\code{\link{find_bicliques}}.}
}
\value{
An object of class \code{"mwperm"}: \code{estimate}/\code{se_naive} are the
OLS estimate and naive SE, \code{conf_int} (or \code{conf_region}/
\code{conf_box} for several coefficients) the IPT inverted-test confidence
set, and \code{pvalue} the IPT permutation p-value; see
\code{\link{mwperm_dyadic}} for the field provenance in full. Its
\code{n_obs} is the number of cells actually used (inside the selected
blocks), and the \code{note} field records how many cells and blocks were
retained.
}
\description{
Finite-sample valid test of \eqn{H_0: \beta = b} in the dyadic model when the
design is \emph{not} a complete array, i.e. some \eqn{(i, j)} cells are
unobserved (Procedure 2 of Guo, Toulis and Wang, 2026). General row/column
permutations are infeasible once cells are missing, so the test is built on
\emph{fully observed bicliques}: maximal-ish blocks of row clusters
\eqn{I_q} and column clusters \eqn{J_q} for which every cell
\eqn{I_q \times J_q} is observed. The block-cyclic permutation group
(\code{\link{build_perm_set}}) is applied within each block, the blocks are
chosen to be disjoint in both rows and columns so that the joint permutation
is a genuine relabelling of exchangeable clusters, and the residual
statistics are pooled across blocks. Cells outside every selected block are
discarded; this is the price of validity under missingness and is reported in
the result.
}
\details{
Finding the largest fully observed biclique is NP-hard, so a greedy
heuristic is used; an approximate (sub-maximal) solution still yields a valid
test, only a less powerful one. The achievable permutation-group order
\code{K + 1} is capped by the smallest selected block, so there is a tension
between including many small blocks (more pooled data, coarser p-value
resolution) and a few large balanced blocks (finer resolution); tune it with
\code{min_block} and \code{K}.
}
\examples{
data(trade_dyadic)
## drop self-trade and a few random links to create missingness
set.seed(1)
d <- trade_dyadic[trade_dyadic$importer != trade_dyadic$exporter, ]
d <- d[sample(nrow(d), round(0.9 * nrow(d))), ]
fit <- with(d, mwperm_missing(y = log_trade, d = log_dist,
                              x = cbind(log_gdp_i, log_gdp_j),
                              row = importer, col = exporter,
                              min_block = 3, seed = 1))
fit
}
\references{
Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference under
multi-way clustering and missing data. arXiv:2601.08610.
}
\seealso{
\code{\link{mwperm_dyadic}} for the complete-array case;
\code{\link{find_bicliques}} for the block finder used internally.
}
