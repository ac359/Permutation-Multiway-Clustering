% Generated manually (no roxygen2/CRAN access in build environment).
\name{build_perm_set}
\alias{build_perm_set}
\title{Construct a random block-cyclic permutation group (Algorithm 1)}
\usage{
build_perm_set(n, K, seed = NULL)
}
\arguments{
\item{n}{Integer, the number of indices to permute (the cluster count along
one dimension).}

\item{K}{Integer, the number of \emph{non-identity} permutations. The group
has order \code{K + 1}, so the smallest attainable p-value is
\code{1 / (K + 1)}. Requires \code{K + 1 <= n}.}

\item{seed}{Optional integer seed for the random relabelling. If \code{NULL},
the current RNG state is used.}
}
\value{
A list of length \code{K + 1} of integer vectors, each a permutation
of \code{seq_len(n)} given as an image vector (entry \code{i} is the image
of \code{i}). The first element is the identity. The list carries the
attribute \code{"block_size"} equal to \code{K + 1}.
}
\description{
Builds a set of \code{K + 1} permutations of \code{seq_len(n)} that form a
cyclic algebraic group, following Algorithm 1 of Guo, Toulis and Wang
(2026). A random one-to-one relabelling
\eqn{\pi} is drawn, the index set is split into consecutive blocks of size
\code{K + 1}, and each non-identity element cyclically shifts every block by
\eqn{k} positions. Because the elements are powers of a single generator,
the returned set is closed under composition (a cyclic group of order
\code{K + 1}).
}
\details{
The randomness in the relabelling is what makes the resulting test a
\emph{random} invariant test: different seeds yield slightly different
p-values. Aggregating several runs (see the \code{n_reps} argument of the
\code{mwperm_*} functions) by taking the median p-value is recommended.
}
\examples{
G <- build_perm_set(n = 8, K = 3, seed = 1)
length(G)          # 4 permutations (identity + 3)
G[[1]]             # identity
attr(G, "block_size")
}
\references{
Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference under
multi-way clustering and missing data. arXiv:2601.08610.
}
