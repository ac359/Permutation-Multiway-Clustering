% Generated manually (no roxygen2/CRAN access in build environment).
\name{confint.mwperm}
\alias{confint.mwperm}
\title{Confidence interval from an inverted permutation test}
\usage{
\method{confint}{mwperm}(object, parm, level = NULL, ...)
}
\arguments{
\item{object}{An object of class \code{"mwperm"}.}

\item{parm}{Ignored; kept for generic compatibility.}

\item{level}{Confidence level; must match the level used at fitting, otherwise
an error is raised (the set cannot be re-derived without the stored
permutations). Defaults to the stored level.}

\item{...}{Ignored.}
}
\value{
A matrix with the lower and upper limits, one row per coefficient: the
inverted-test interval for a single coefficient, or the marginal extent of the
joint confidence region for several (see \code{object$conf_region} for the full
set of retained vectors).
}
\description{
Extracts the test-inversion confidence set stored in a \code{"mwperm"} object:
the interval for a single coefficient, or the marginal extent of the joint
confidence region for several. Note that \code{level} is fixed at fitting time
(\code{1 - alpha}); a different level requires refitting with the corresponding
\code{alpha}.
}
\seealso{
\code{\link{mwperm_dyadic}}, \code{\link{print.mwperm}}.
}
