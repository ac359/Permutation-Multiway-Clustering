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
\code{conf_high} and \code{p_value}. \code{estimate} is the \emph{OLS} point
estimate (least squares on the full design); \code{conf_low}/\code{conf_high}
are the limits of the IPT (inverted permutation test) confidence set -- the
inverted-test interval for a single coefficient, or the marginal extent of
the joint confidence region for several; \code{p_value} is the permutation
p-value of the joint test (repeated across rows when there are several
coefficients).
}
\description{
Returns (invisibly, after printing) a one-row-per-coefficient data frame with
the OLS point estimate, naive standard error, the IPT confidence limits, and
the permutation p-value for \eqn{H_0:\beta = b}; see \strong{Value} for the
exact meaning of each column.
}
\seealso{
\code{\link{print.mwperm}}, \code{\link{confint.mwperm}}.
}
