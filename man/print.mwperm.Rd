% Generated manually (no roxygen2/CRAN access in build environment).
\name{print.mwperm}
\alias{print.mwperm}
\title{Print a multi-way permutation test}
\usage{
\method{print}{mwperm}(x, digits = 4, ...)
}
\arguments{
\item{x}{An object of class \code{"mwperm"}.}

\item{digits}{Number of significant digits for the estimate and interval.}

\item{...}{Ignored.}
}
\value{
\code{x}, invisibly.
}
\description{
Compact display of the design, the permutation-group order, the point
estimate and inverted confidence interval (when available), and the
permutation p-value with the reject/retain decision.
}
\seealso{
\code{\link{mwperm_dyadic}}, \code{\link{summary.mwperm}}.
}
