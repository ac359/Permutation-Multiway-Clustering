% Generated manually (no roxygen2/CRAN access in build environment).
\name{plot.mwperm}
\alias{plot.mwperm}
\title{Visualise the permutation reference distribution}
\usage{
\method{plot}{mwperm}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"mwperm"}.}

\item{...}{Passed to the underlying plotting calls.}
}
\value{
\code{x}, invisibly.
}
\description{
Plots the per-replication p-values when several replications were run, with
the reported (median) p-value and the level \code{alpha} marked; for a single
replication it draws a simple bar of the p-value against \code{alpha}. This
gives a quick read on the Monte-Carlo stability of the randomised test.
}
\seealso{
\code{\link{mwperm_dyadic}}, \code{\link{print.mwperm}}.
}
