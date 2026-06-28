% Generated manually (no roxygen2/CRAN access in build environment).
\name{summary.mwperm}
\alias{summary.mwperm}
\title{Summarise a multi-way permutation test}
\usage{
\method{summary}{mwperm}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"mwperm"}.}

\item{...}{Ignored.}
}
\value{
A data frame, invisibly, with one row per coefficient and columns
\code{term}, \code{estimate}, \code{se_naive}, \code{conf_low},
\code{conf_high} and \code{p_value}.
}
\description{
Returns (invisibly, after printing) a one-row-per-coefficient data frame with
the point estimate, naive standard error, and -- for the first coefficient --
the inverted confidence-interval end points, alongside the permutation
p-value for \eqn{H_0:\beta = b}.
}
\seealso{
\code{\link{print.mwperm}}, \code{\link{confint.mwperm}}.
}
