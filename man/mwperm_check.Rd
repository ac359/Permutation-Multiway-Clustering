% Generated manually (no roxygen2/CRAN access in build environment).
\name{mwperm_check}
\alias{mwperm_check}
\alias{print.mwperm_design}
\title{Diagnose a multi-way clustered dataset and choose the appropriate design}
\usage{
mwperm_check(index, y = NULL, d = NULL, data = NULL,
             time = NULL, rep = NULL,
             design = c("auto", "dyadic", "threeway", "panel",
                        "layout", "missing"))

\method{print}{mwperm_design}(x, ...)
}
\arguments{
\item{index}{The clustering dimensions (2 or 3): a data frame, a named list
of vectors, or a character vector of column names resolved against
\code{data}.}

\item{y, d}{Optional outcome and covariate(s) of interest; only used for
extra diagnostics (e.g. the layout no-power warning when \code{d} is
constant within every cell), never for fitting.}

\item{data}{Optional data frame against which character \code{index},
\code{time} and \code{rep} entries are resolved.}

\item{time}{Optional explicit time dimension: a vector, or the name of a
column of \code{data} (or of one of the \code{index} columns). Forces the
panel interpretation of that dimension.}

\item{rep}{Optional explicit replication identifier (vector or column name):
declares within-cell replication and forces the layout design.}

\item{design}{Force a design instead of auto-detecting (the structure is
still validated against it).}

\item{x}{An object of class \code{"mwperm_design"} (print method).}

\item{...}{Ignored.}
}
\value{
An object of class \code{"mwperm_design"}: a list with fields \code{design}
(the chosen design), \code{roles} (which index plays row/col/id1..3/time/rep),
\code{dims} (levels per dimension), \code{n_obs}, \code{cells}
(observed/expected), \code{balance}, \code{K_default} and
\code{resolution_ok} (whether a 95\% confidence set is attainable),
\code{call_str} (the downstream call), \code{reason} (one-line explanation),
and \code{warnings}/\code{notes} (the assumption-fork notices etc.). Its
\code{print} method lays this out as a short human diagnosis.
}
\description{
Inspects the clustering structure of a dataset -- number of index dimensions,
repeated cells, completeness/balance -- and reports which \code{mwperm_*}
test applies, without running any permutations or fitting anything.
\code{\link{mwperm}} uses it for automatic dispatch; call it directly to see
the diagnosis.
}
\details{
Structural forks (dyadic vs missing vs layout-by-replication) are resolved
silently from the data. Two forks depend on an \emph{exchangeability
assumption the data cannot reveal} and are therefore announced with override
instructions, defaulting to the choice that remains valid under the widest
set of error processes:
\itemize{
  \item \strong{panel vs three-way} (complete balanced 3-index arrays):
    running \code{\link{mwperm_threeway}} on a panel whose errors are
    dependent over time is \emph{invalid} (size distortion), while running
    \code{\link{mwperm_panel}} on genuinely three-way exchangeable data is
    valid, only less powerful. The default is therefore \code{panel}; a
    time-like third index (by name, temporal class, or regular numeric
    spacing) is used as the time dimension silently, an ambiguous one with a
    notice. Force \code{design = "threeway"} only when all three dimensions
    are genuinely exchangeable.
  \item \strong{layout vs suppressed panel} (2 indices with repeated cells):
    repeats are treated as within-cell replication
    (\code{\link{mwperm_layout}}), which assumes the replicates are
    exchangeable within cells -- if they are really a time series, pass the
    time variable via \code{time =} to get the panel test instead. A notice
    is attached.
}
}
\examples{
data(trade_dyadic)
mwperm_check(index = c("importer", "exporter"), data = trade_dyadic)
data(trade_panel)
mwperm_check(index = c("importer", "exporter", "year"), data = trade_panel)
}
\seealso{
\code{\link{mwperm}} for one-call dispatch.
}
