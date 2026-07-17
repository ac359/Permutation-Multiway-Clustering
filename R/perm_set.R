#' Construct a random block-cyclic permutation group (Algorithm 1)
#'
#' Builds a set of \code{K + 1} permutations of \code{seq_len(n)} that form a
#' cyclic algebraic group, following Algorithm 1 of Guo, Toulis and Wang
#' (2026). A random one-to-one relabelling
#' \eqn{\pi} is drawn, the index set is split into consecutive blocks of size
#' \code{K + 1}, and each non-identity element cyclically shifts every block by
#' \eqn{k} positions. Because the elements are powers of a single generator,
#' the returned set is closed under composition (a cyclic group of order
#' \code{K + 1}).
#'
#' The randomness in the relabelling is what makes the resulting test a
#' \emph{random} invariant test: different seeds yield slightly different
#' p-values. Aggregating several runs (see the \code{n_reps} argument of the
#' \code{mwperm_*} functions) by taking the median p-value is recommended;
#' see \code{\link{mwperm_dyadic}} for the exact status of that aggregation.
#'
#' \strong{Reproducibility.} A seeded result is reproducible only under the
#' same RNG configuration: the same \code{RNGkind()} (generator \emph{and}
#' sample kind -- R's defaults changed in 3.6.0) and, when cluster ids are
#' supplied as character strings, the same collation locale
#' (\code{LC_COLLATE} determines factor level order and hence the dense id
#' coding the relabelling acts on). Integer or factor cluster ids make seeded
#' results locale-proof. Validity is unaffected either way: whatever group is
#' realised is a genuine cyclic group, so the test is exact under any RNG
#' configuration -- only cross-environment reproducibility depends on it.
#'
#' \strong{What certifies the group property.} Closure under composition --
#' the property Theorem 1 of the paper rests on -- is verified algebraically
#' by the shipped test suite (full composition tables, including at the
#' observation level for every design). Monte-Carlo size simulations cannot
#' detect closure defects, because a broken set whose elements are still
#' per-dimension permutations continues to control size element-wise;
#' simulation evidence therefore never certifies the group structure.
#'
#' @param n Integer, the number of indices to permute (the cluster count along
#'   one dimension).
#' @param K Integer, the number of \emph{non-identity} permutations. The group
#'   has order \code{K + 1}, so the smallest attainable p-value is
#'   \code{1 / (K + 1)}. Requires \code{K + 1 <= n}.
#' @param seed Optional integer seed for the random relabelling. If \code{NULL},
#'   the current RNG state is used.
#'
#' @return A list of length \code{K + 1} of integer vectors, each a permutation
#'   of \code{seq_len(n)} given as an image vector (entry \code{i} is the image
#'   of \code{i}). The first element is the identity. The list carries the
#'   attribute \code{"block_size"} equal to \code{K + 1}.
#'
#' @references
#' Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference under
#' multi-way clustering and missing data. arXiv:2601.08610.
#'
#' @examples
#' G <- build_perm_set(n = 8, K = 3, seed = 1)
#' length(G)          # 4 permutations (identity + 3)
#' G[[1]]             # identity
#' attr(G, "block_size")
#' @export
build_perm_set <- function(n, K, seed = NULL) {
  n <- as.integer(n)
  K <- as.integer(K)
  if (n < 2L) stop("`n` must be at least 2.", call. = FALSE)
  if (K < 1L) stop("`K` must be at least 1.", call. = FALSE)
  B <- K + 1L                      # block size = group order
  if (B > n) {
    stop(sprintf(
      "K + 1 = %d exceeds n = %d: no non-trivial permutation exists. Use K <= n - 1.",
      B, n), call. = FALSE)
  }

  ## Reproducible relabelling without disturbing the caller's RNG stream: save
  ## the global seed, reseed, and restore on exit (see .save_seed/.restore_seed).
  if (!is.null(seed)) {
    old <- .save_seed()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }

  ## Random relabelling pi (a permutation of 1..n) and its inverse. pi is what
  ## makes the test a *random* invariant test; different seeds give different pi.
  pi_vec <- sample.int(n)           # pi_vec[i] = pi(i)
  pi_inv <- integer(n)              # pi_inv[pi(i)] = i  (the inverse map)
  pi_inv[pi_vec] <- seq_len(n)

  ## Split 1..n into consecutive blocks of size B; any tail of < B leftover
  ## indices is held fixed (it cannot be cyclically shifted within a full block).
  nb <- n %/% B                     # number of full blocks
  in_block <- seq_len(nb * B)       # the indices that live in a full block
  pos0 <- (in_block - 1L) %% B      # 0-based position within the block
  blk_start <- in_block - pos0      # 1-based index of the block's first element

  ## Each non-identity element k cyclically shifts every block by k positions;
  ## element 0 is the identity. Building all powers of one generator guarantees
  ## the set is a cyclic group of order B (closed under composition).
  perms <- vector("list", B)
  for (k in 0:K) {
    psit <- seq_len(n)              # psi-tilde_k as an image vector; default = identity
    if (k > 0L && nb > 0L) {
      new_pos0 <- (pos0 + k) %% B   # shifted position within the block
      psit[in_block] <- blk_start + new_pos0
    }
    ## Conjugate the block shift by the relabelling: psi_k(i) = pi^{-1}( psi-tilde_k( pi(i) ) )
    perms[[k + 1L]] <- pi_inv[psit[pi_vec]]
  }
  attr(perms, "block_size") <- B
  perms
}

## ---- internal RNG-state helpers ------------------------------------------
## These let a seeded routine draw random numbers reproducibly without
## clobbering the caller's random stream: snapshot .Random.seed, reseed, do the
## work, then restore the snapshot on exit.

#' Snapshot the current global RNG state (or NULL if none has been initialised).
#' @keywords internal
#' @noRd
.save_seed <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
}

#' Restore a global RNG state previously captured by \code{.save_seed()}.
#' A NULL snapshot means the RNG was uninitialised, so we remove the seed again.
#' @keywords internal
#' @noRd
.restore_seed <- function(old) {
  if (is.null(old)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  } else {
    assign(".Random.seed", old, envir = .GlobalEnv)
  }
  invisible(NULL)
}
