% Generated manually (no roxygen2/CRAN access in build environment).
\name{mwperm-package}
\alias{mwperm-package}
\alias{mwperm}
\docType{package}
\title{mwperm: Invariant Permutation Tests for Multi-Way Clustered and Panel Regression}
\description{
Finite-sample valid tests and confidence intervals for regression
coefficients under multi-way (e.g. dyadic) clustering, including panel and
missing-data designs, implementing the invariant permutation test of Guo,
Toulis and Wang (2026) building on Wen, Wang and Wang (2025).
}
\details{
The key object is the block-cyclic permutation group built by
\code{\link{build_perm_set}}. The user-facing tests are
\code{\link{mwperm_dyadic}} (two-way / dyadic clustering),
\code{\link{mwperm_threeway}} (three-way clustering),
\code{\link{mwperm_panel}} (panels with an arbitrary time effect),
\code{\link{mwperm_layout}} (replicated two-way layouts) and
\code{\link{mwperm_missing}} (incomplete arrays, via fully observed
bicliques). All return an object of class \code{"mwperm"} with
\code{\link{print.mwperm}}, \code{\link{summary.mwperm}},
\code{\link{confint.mwperm}} and \code{\link{plot.mwperm}} methods.

Two synthetic data sets, \code{\link{trade_dyadic}} and
\code{\link{trade_panel}}, illustrate the dyadic and panel work flows.
}
\references{
Guo, F. R., Toulis, P. and Wang, Y. (2026). Permutation inference under
multi-way clustering and missing data.

Wen, K., Wang, T. and Wang, Y. (2025). Residual permutation test for
regression coefficient testing. \emph{The Annals of Statistics} 53(2),
724--748.
}
\seealso{
\code{\link{mwperm_dyadic}}, \code{\link{mwperm_panel}},
\code{\link{build_perm_set}}.
}
\keyword{package}
