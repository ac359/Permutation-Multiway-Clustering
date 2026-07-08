% Generated manually (no roxygen2/CRAN access in build environment).
\name{mwperm_threeway}
\alias{mwperm_threeway}
\title{Invariant permutation test under three-way clustering}
\usage{
mwperm_threeway(y, d, x = NULL, id1, id2, id3, K = NULL,
                alpha = 0.05, beta_null = 0, conf_int = TRUE,
                n_reps = 1L, seed = NULL, grid = NULL, n_cores = 1L)
}
\arguments{
\item{y}{Numeric outcome vector, one entry per observed cell.}

\item{d}{Numeric vector or matrix of the covariate(s) of interest. With a
single covariate a confidence interval is produced; with several, a joint
confidence region.}

\item{x}{Optional numeric matrix or data frame of nuisance covariates; an
intercept is always added internally. May be \code{NULL}.}

\item{id1, id2, id3}{Cluster identifiers for the three dimensions.}

\item{K}{Number of non-identity permutations; defaults to
\code{min(m, n, ell) - 1} capped at 199.}

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
Finite-sample valid test of \eqn{H_0: \beta = b} in the three-way model
\deqn{y_{ijl} = x_{ijl}^\top \gamma + d_{ijl}^\top \beta + \varepsilon_{ijl},}
\eqn{i \in [m], j \in [n], l \in [\ell]}, under full three-way
exchangeability of the errors (condition InvA of Guo, Toulis and Wang,
2026),
\deqn{(\varepsilon_{ijl}) \;{\buildrel d \over =}\;
      (\varepsilon_{\pi(i)\sigma(j)\psi(l)}) \mid X, D,}
which holds, e.g., for the random-effects structure \eqn{\varepsilon_{ijl} =
\eta_i + \xi_j + \zeta_l + u_{ijl}} with all components i.i.d. within family.
A permutation is drawn independently for each of the three dimensions and
applied jointly.
}
\details{
Use this when all three dimensions are genuinely exchangeable (e.g.
importer, exporter and product category). If one dimension is time with
autocorrelation, use \code{\link{mwperm_panel}} instead. The data must form
a complete balanced array.
}
\examples{
## small balanced 3-way array with a random-effects error
set.seed(1)
m <- 8
g <- expand.grid(i = seq_len(m), j = seq_len(m), l = seq_len(m))
a <- rnorm(m); b <- rnorm(m); cc <- rnorm(m)
g$d <- rnorm(nrow(g))
g$y <- 0.6 * g$d + a[g$i] + b[g$j] + cc[g$l] + rnorm(nrow(g))
fit <- with(g, mwperm_threeway(y = y, d = d, id1 = i, id2 = j, id3 = l,
                               conf_int = FALSE, seed = 1))
fit
}
\references{
Guo, F. R., Toulis, P. and Wang, Y. (2026). Permutation inference under
multi-way clustering and missing data, Section 6.1.
}
\seealso{
\code{\link{mwperm_dyadic}}, \code{\link{mwperm_panel}}.
}
