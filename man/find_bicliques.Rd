% Generated manually (no roxygen2/CRAN access in build environment).
\name{find_bicliques}
\alias{find_bicliques}
\title{Disjoint fully observed bicliques (greedy or exact)}
\usage{
find_bicliques(row, col, min_block = 2L,
               method = c("greedy", "exact"), node_budget = 2e5L)
}
\arguments{
\item{row, col}{Integer (or factor-coercible) cluster ids of the observed
cells; the two vectors have equal length, one entry per observed cell.}

\item{min_block}{Integer; blocks whose smaller side is below this are not
returned. Defaults to 2.}

\item{method}{Either \code{"greedy"} (default) or \code{"exact"}; see Details.}

\item{node_budget}{Integer node cap for the \code{"exact"} branch-and-bound
per block (default 200000). Ignored when \code{method = "greedy"}.}
}
\value{
A list of blocks, each a list with integer components \code{rows} and
\code{cols} giving the (original-coding) cluster ids in that biclique.
}
\description{
Partitions (part of) a sparse two-way layout into blocks of row clusters
\eqn{I_q} and column clusters \eqn{J_q} such that every cell of
\eqn{I_q \times J_q} is observed and the blocks share no row or column
cluster. Blocks are peeled off one at a time: the largest fully observed
biclique among the still-available rows and columns is found, removed, and
the search repeats on the remainder. Used by \code{\link{mwperm_missing}}.
}
\details{
Two ways to find each block are offered through \code{method}:
\describe{
  \item{\code{"greedy"}}{(default) a fast seed-and-intersect heuristic that,
    starting from each of a few highest-degree rows, adds rows while keeping
    the common set of fully observed columns. Always valid but may return a
    sub-maximal block.}
  \item{\code{"exact"}}{a branch-and-bound search for the maximum-area fully
    observed biclique, with area-based pruning. Maximum-edge biclique is
    NP-hard, so a per-block node budget caps the work; if it is hit the search
    falls back to the greedy block for that step (still valid) and a warning
    is issued. For the modest cluster counts this method targets the exact
    search is typically fast and returns blocks at least as large as greedy,
    giving finer p-value resolution.}
}
}
\examples{
## a 4x4 grid missing its diagonal
g <- expand.grid(i = 1:4, j = 1:4)
g <- g[g$i != g$j, ]
find_bicliques(g$i, g$j, min_block = 2)
find_bicliques(g$i, g$j, min_block = 2, method = "exact")
}
\seealso{
\code{\link{mwperm_missing}}.
}
