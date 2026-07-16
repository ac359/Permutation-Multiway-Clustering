% Generated manually (no roxygen2/CRAN access in build environment).
\name{mwperm_panel}
\alias{mwperm_panel}
\title{Invariant permutation test for panel (longitudinal) dyadic regression}
\usage{
mwperm_panel(y, d, x = NULL, row, col, time, K = NULL,
             alpha = 0.05, beta_null = 0, conf_int = TRUE,
             n_reps = 1L, seed = NULL, grid = NULL, time_fe = TRUE,
             n_cores = 1L)
}
\arguments{
\item{y}{Numeric outcome vector, one entry per observed cell.}

\item{d}{Numeric vector or matrix of the covariate(s) of interest
\eqn{d_{ijt}} (may be time-varying). With a single covariate a confidence
interval is produced; with several, a joint confidence region.}

\item{x}{Optional nuisance covariates \eqn{x_{ijt}}; intercept added
internally. May be \code{NULL}.}

\item{row, col}{Row- and column-cluster identifiers (length = number of
observations).}

\item{time}{Time-period identifiers (length = number of observations).}

\item{K}{Number of non-identity permutations; defaults to
\code{min(n_row, n_col) - 1} capped at 199.}

\item{alpha}{Significance level for the decision and the inverted interval.}

\item{beta_null}{Null value \eqn{b} to test.}

\item{conf_int}{Logical; compute a test-inversion confidence set: an interval
for a single covariate, or a joint confidence region for several.}

\item{n_reps}{Number of independent runs aggregated by the median.}

\item{seed}{Optional integer; run \code{r} uses seed \code{seed + r - 1}.}

\item{grid}{Optional numeric vector of candidate \eqn{\beta} values for the
confidence interval.}
\item{n_cores}{Number of CPU cores for the permutation computations
(default 1 = serial). Parallelism is over the \code{n_reps} repetitions when
several are run with a \code{seed}, otherwise over the \code{K}
per-permutation factorizations; either way the result is \emph{identical} to
the serial one. See \code{\link{mwperm_dyadic}} for details.}

\item{time_fe}{Logical; if \code{TRUE} (the default) time fixed effects are
added to the nuisance design. Under condition InvB this is valid (the time
dummies are invariant to the within-period permutation) and it removes the
time trend \eqn{\zeta_t} from the residuals, which de-biases the reported
point estimate and sharpens the test when treatment timing is correlated
with the period.}
}
\value{
An object of class \code{"mwperm"}: \code{estimate}/\code{se_naive} are the
OLS estimate and naive SE, \code{conf_int} (or \code{conf_region}/
\code{conf_box} for several coefficients) the IPT inverted-test confidence
set, and \code{pvalue} the IPT permutation p-value; see
\code{\link{mwperm_dyadic}} for the field provenance in full.
}
\description{
Finite-sample valid test of \eqn{H_0: \beta = b} in the panel model
\deqn{y_{ijt} = x_{ijt}^\top \gamma + d_{ijt}^\top \beta + \varepsilon_{ijt},}
where \eqn{i, j} index two cross-sectional clustering dimensions (e.g.
importer and exporter countries) and \eqn{t} indexes time. Full three-way
exchangeability is implausible because errors are typically autocorrelated
over time. Instead this test only assumes exchangeability across the first
two dimensions \emph{within} each time period (condition InvB of Guo, Toulis
and Wang, 2026),
\deqn{(\varepsilon_{ijt})_{i,j} \;{\buildrel d \over =}\;
      (\varepsilon_{\pi(i)\sigma(j)t})_{i,j} \mid X, D,}
which holds, for instance, under \eqn{\varepsilon_{ijt} = \eta_i + \xi_j +
\zeta_t + u_{ijt}} with an \emph{arbitrary} time trend \eqn{\zeta_t}. The
same row/column permutation is applied in every period, so any unknown time
effect is held fixed and partialled out. This is, to the authors' knowledge,
the first finite-sample valid test for \eqn{\beta = 0} under exchangeable
errors in such panel models.
}
\details{
The data must form a complete balanced array: every \eqn{(i, j, t)} cell
present exactly once.
}
\examples{
data(trade_panel)
fit <- with(trade_panel,
            mwperm_panel(y = log_trade, d = fta, x = cbind(log_gdp_i, log_gdp_j),
                         row = importer, col = exporter, time = year, seed = 1))
fit
}
\references{
Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference under
multi-way clustering and missing data, Section 6.2. arXiv:2601.08610.
}
\seealso{
\code{\link{mwperm_dyadic}}, \code{\link{mwperm_threeway}}.
}
