% Generated manually (no roxygen2/CRAN access in build environment).
\name{mwperm}
\alias{mwperm}
\title{One-call invariant permutation test with automatic design detection}
\usage{
mwperm(y, d, x = NULL, index, data = NULL, time = NULL, rep = NULL,
       design = c("auto", "dyadic", "threeway", "panel", "layout", "missing"),
       K = NULL, alpha = 0.05, beta_null = 0, conf_int = TRUE,
       n_reps = 1L, seed = NULL, grid = NULL, n_cores = 1L,
       time_fe = TRUE, L0 = NULL, min_block = 3L,
       block_method = c("greedy", "exact"), verbose = TRUE)
}
\arguments{
\item{y, d, x}{Outcome, covariate(s) of interest, and optional nuisance
covariates, as in \code{\link{mwperm_dyadic}}. With \code{data} given, each
may also be a character (vector of) column name(s) resolved against it.}

\item{index}{The clustering dimensions (2 or 3): a data frame, a named list
of vectors, or a character vector of column names in \code{data}.}

\item{data}{Optional data frame; column names in \code{y}, \code{d},
\code{x}, \code{index}, \code{time}, \code{rep} are resolved against it.}

\item{time, rep}{Optional explicit role tags (vector or column name); see
\code{\link{mwperm_check}}.}

\item{design}{Force a design instead of auto-detecting (the structure is
still validated against it).}

\item{K}{Number of non-identity permutations; see the dispatched function.}

\item{alpha}{Significance level (\code{conf_level = 1 - alpha}).}

\item{beta_null}{Null value \eqn{b} to test.}

\item{conf_int}{Logical; compute a confidence set by test inversion.}

\item{n_reps}{Number of independent runs aggregated by the median.}

\item{seed}{Optional integer seed; run \code{r} uses \code{seed + r - 1}.}

\item{grid}{Optional candidate \eqn{\beta} values for the confidence set.}

\item{n_cores}{Number of CPU cores; see \code{\link{mwperm_dyadic}}.}

\item{time_fe}{Passed to \code{\link{mwperm_panel}} (panel only; supplying
it for another design warns and ignores it).}

\item{L0}{Passed to \code{\link{mwperm_layout}} (layout only).}

\item{min_block, block_method}{Passed to \code{\link{mwperm_missing}}
(missing only).}

\item{verbose}{If \code{TRUE} (default) print one line stating the detected
design and the dispatched call.}
}
\value{
The \code{"mwperm"} object of the dispatched test, with an extra \code{auto}
field recording the detection (design, roles, reason); the detection notices
are prepended to the object's \code{note} field and shown by
\code{\link{print.mwperm}}, and any assumption-fork or weak-evidence notice
is additionally raised as a \code{warning} at fit time. Field provenance
(see \code{\link{mwperm_dyadic}}
for the full account): \code{estimate}/\code{se_naive} are the OLS estimate
and naive SE, \code{conf_int} (or \code{conf_region}/\code{conf_box} for
several coefficients) the IPT inverted-test confidence set, and
\code{pvalue} the IPT permutation p-value.
}
\description{
Detects the clustering design of the data via \code{\link{mwperm_check}} and
dispatches to the matching test -- \code{\link{mwperm_dyadic}},
\code{\link{mwperm_panel}}, \code{\link{mwperm_threeway}},
\code{\link{mwperm_layout}} or \code{\link{mwperm_missing}} -- forwarding all
arguments unchanged. A thin convenience layer: the returned object is exactly
what the underlying function returns (plus a record of what was detected),
and calling the specific function directly with the same seed gives identical
results.
}
\details{
See \code{\link{mwperm_check}} for the detection rules, in particular the two
assumption-dependent forks (panel-vs-threeway and layout-vs-suppressed-panel)
that are announced rather than silently resolved. Structural forks (complete
vs incomplete arrays, replicated cells) are resolved silently.
}
\examples{
data(trade_dyadic)
fit <- mwperm(y = "log_trade", d = "log_dist",
              x = c("log_gdp_i", "log_gdp_j"),
              index = c("importer", "exporter"),
              data = trade_dyadic, seed = 1)
fit
}
\seealso{
\code{\link{mwperm_check}} for the diagnosis without any computation;
\code{\link{mwperm_dyadic}}, \code{\link{mwperm_panel}},
\code{\link{mwperm_threeway}}, \code{\link{mwperm_layout}},
\code{\link{mwperm_missing}} for the underlying tests.
}
