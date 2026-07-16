% Generated manually (no roxygen2/CRAN access in build environment).
\name{mwperm_layout}
\alias{mwperm_layout}
\title{Invariant permutation test for two-way layouts (within-cell replication)}
\usage{
mwperm_layout(y, d, x = NULL, row, col, rep = NULL, L0 = NULL, K = NULL,
              alpha = 0.05, beta_null = 0, conf_int = TRUE,
              n_reps = 1L, seed = NULL, grid = NULL, n_cores = 1L)
}
\arguments{
\item{y}{Numeric outcome vector, one entry per observation.}

\item{d}{Numeric vector or matrix of the covariate(s) of interest. Must vary
within cells (see Details). With a single covariate a confidence interval is
produced; with several, a joint confidence region.}

\item{x}{Optional numeric matrix or data frame of nuisance covariates; an
intercept is always added internally. May be \code{NULL}.}

\item{row, col}{Cell identifiers along the two layout dimensions.}

\item{rep}{Optional replication identifier within each cell; if \code{NULL},
the order of appearance within a cell is used.}

\item{L0}{Optional integer capacity threshold for balancing an unbalanced
layout (Section 6.3). When supplied, cells with fewer than \code{L0}
replicates are dropped and each remaining cell is uniformly downsampled to
exactly \code{L0} replicates, giving a balanced array before testing. The
downsampling is reproducible through \code{seed}. When \code{NULL} (the
default) all replicates are used with their (possibly unequal) cell sizes.}

\item{K}{Number of non-identity permutations; defaults to
\code{min(cell size) - 1} capped at 199 (with \code{L0}, to \code{L0 - 1}).
Must satisfy \code{K + 1 <= min(cell size)}.}

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
}
\value{
An object of class \code{"mwperm"}: \code{estimate}/\code{se_naive} are the
OLS estimate and naive SE, \code{conf_int} (or \code{conf_region}/
\code{conf_box} for several coefficients) the IPT inverted-test confidence
set, and \code{pvalue} the IPT permutation p-value; see
\code{\link{mwperm_dyadic}} for the field provenance in full.
}
\description{
Finite-sample valid test of \eqn{H_0: \beta = b} for a two-way layout
\deqn{y_{ijl} = x_{ijl}^\top \gamma + d_{ijl}^\top \beta + \varepsilon_{ijl}}
in which each \eqn{(i, j)} cell contains \eqn{\ell_{ij}} replicate
observations and the cell sizes may be \emph{unequal}. Permutations are
applied only within cells, over the replication index \eqn{l}. This is valid
under the relaxed condition that errors are exchangeable with respect to
\eqn{l} within each cell (condition InvA restricted to \eqn{l}), e.g.
\eqn{\varepsilon_{ijl} = \eta_{ij} + \zeta_l + u_{ijl}} with arbitrary
\eqn{\eta_{ij}} but i.i.d. \eqn{\zeta_l}. It is appropriate when \eqn{l}
indexes independent replications (two-way layouts in randomised
experiments). See Guo, Toulis and Wang (2026), Section 6.3.
}
\details{
Note: the covariate of interest must vary within cells. If \eqn{d_{ijl}} is
constant within every \eqn{(i, j)} cell (a cell-level covariate), within-cell
permutation yields a trivial test with no power; a warning is issued.

Supplying \code{L0} balances an unbalanced layout to exactly \code{L0}
replicates per cell (cells with fewer are dropped, denser cells are
uniformly downsampled). A uniform random subset of exchangeable replicates is
itself exchangeable, so validity is preserved while the permutation-group
order becomes the same (\code{L0}) in every cell.
}
\examples{
## 6 x 6 cells, 6 within-cell replicates, treatment varies within cell
set.seed(1)
cells <- expand.grid(i = 1:6, j = 1:6)
dat <- do.call(rbind, lapply(seq_len(nrow(cells)), function(r) {
  data.frame(i = cells$i[r], j = cells$j[r], l = 1:6,
             d = rnorm(6),
             eta = rnorm(1))
}))
dat$y <- 0.6 * dat$d + dat$eta + rnorm(nrow(dat))
fit <- with(dat, mwperm_layout(y = y, d = d, row = i, col = j, rep = l,
                               n_reps = 5, conf_int = FALSE, seed = 1))
fit
}
\references{
Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference under
multi-way clustering and missing data, Section 6.3. arXiv:2601.08610.
}
\seealso{
\code{\link{mwperm_dyadic}}, \code{\link{mwperm_missing}}.
}
