## ============================================================================
## R/perm_set.R -- Algorithm 1: the random block-cyclic permutation group
##
## Purpose. build_perm_set() draws, for ONE clustering dimension, the K + 1
##   permutations {psi_0 = Id, psi_1, ..., psi_K} that Procedure 1 of Guo,
##   Toulis & Wang (2026) ("GTW") compares the data against; .save_seed() and
##   .restore_seed() keep every seeded draw in the package off the caller's
##   random stream.
## Paper. GTW Algorithm 1 (Section 3.2) and Proposition 2 (closure under
##   composition, which Theorem 1 needs). Adapted from Algorithm 1 of Wen,
##   Wang & Wang (2025, "WWW"); like GTW, it omits WWW's repeat-until loop on
##   the trace criterion, which serves power, not validity.
## Two readings of the printed algorithm, both recorded in TESTING_PLAN.md:
##   * the display's "i mod (K+1)" must be the 1-based residue
##     ((i - 1) mod (K + 1)) + 1, or psi~_k is not a permutation (D1);
##   * the worked example S_1^k lists the inverse of the displayed map (D2).
##   The code follows the display, so element k + 1 of the returned list is
##   psi_k; the example's ordering would give the same SET.
## Pipeline. mwperm() -> dispatch -> design worker -> [permutation
##   construction] -> projection engine -> median aggregation -> test
##   inversion -> S3 methods. Every permutation front end's perm_builder
##   calls build_perm_set() once per permuted dimension (or cell, or block)
##   per repetition, and the gather-vector builders in core.R, layout.R and
##   missing.R lift the result to the observations.
## ============================================================================

#' Construct a random block-cyclic permutation group (Algorithm 1)
#'
#' Builds a set of `K + 1` permutations of `seq_len(n)` that form a cyclic
#' algebraic group, following Algorithm 1 of Guo, Toulis and Wang (2026). A
#' random one-to-one relabelling pi is drawn, the index set is split into
#' consecutive blocks of size `K + 1`, and each non-identity element
#' cyclically shifts every block by k positions. Because the elements are
#' powers of a single generator, the returned set is closed under composition
#' (a cyclic group of order `K + 1`).
#'
#' The randomness in the relabelling is what makes the resulting test a
#' *random* invariant test: different seeds yield slightly different p-values.
#' Aggregating several runs (see the `n_reps` argument of the `mwperm_*`
#' functions) by taking the median p-value is recommended; see
#' [mwperm_dyadic()] for the exact status of that aggregation.
#'
#' **Reproducibility.** A seeded result is reproducible only under the same
#' RNG configuration: the same `RNGkind()` (generator *and* sample kind -- R's
#' defaults changed in 3.6.0) and, when cluster ids are supplied as character
#' strings, the same collation locale (`LC_COLLATE` determines factor level
#' order and hence the dense id coding the relabelling acts on). Integer or
#' factor cluster ids make seeded results locale-proof. Validity is unaffected
#' either way: whatever group is realised is a genuine cyclic group, so the
#' test is exact under any RNG configuration -- only cross-environment
#' reproducibility depends on it.
#'
#' **What certifies the group property.** Closure under composition -- the
#' property Theorem 1 of the paper rests on -- is verified algebraically by
#' the shipped test suite (full composition tables, including at the
#' observation level for every design). Monte-Carlo size simulations cannot
#' detect closure defects, because a broken set whose elements are still
#' per-dimension permutations continues to control size element-wise;
#' simulation evidence therefore never certifies the group structure.
#'
#' @details Implements Algorithm 1 of Guo, Toulis and Wang (2026). The residue
#'   in its display is read as 1-based, ((i - 1) mod (K + 1)) + 1, the reading
#'   under which psi~_k is a permutation; element k + 1 of the list is
#'   `psi_k = pi^{-1} o psi~_k o pi`, so the list is closed under composition
#'   (their Proposition 2). The paper's worked example lists the same set in
#'   the opposite order: its element k is element K + 1 - k here.
#' @param n Integer, the number of indices to permute (the cluster count along
#'   one dimension).
#' @param K Integer, the number of *non-identity* permutations. The group has
#'   order `K + 1`, so the smallest attainable p-value is `1 / (K + 1)`.
#'   Requires `K + 1 <= n`.
#' @param seed Optional integer seed for the random relabelling. If `NULL`,
#'   the current RNG state is used.
#'
#' @return A list of length `K + 1` of integer vectors, each a permutation of
#'   `seq_len(n)` given as an image vector (entry `i` is the image of `i`).
#'   The first element is the identity. The list carries the attribute
#'   `"block_size"` equal to `K + 1`.
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data. arXiv:2601.08610.
#'
#' @seealso [mwperm_dyadic()] and the other front ends, which compose these
#'   groups into observation-level permutations; [build_flip_set()] for the
#'   package's other group construction, the sign-flip group behind
#'   [mwperm_dyadic_het()].
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
      paste0("K + 1 = %d exceeds n = %d: no non-trivial permutation ",
             "exists. Use K <= n - 1."),
      B, n), call. = FALSE)
  }

  ## Reproducible relabelling without disturbing the caller's RNG stream: save
  ## the global seed, reseed, and restore on exit (see
  ## .save_seed/.restore_seed).
  if (!is.null(seed)) {
    old <- .save_seed()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }

  ## Random relabelling pi (a permutation of 1..n) and its inverse. pi is what
  ## makes the test a *random* invariant test; different seeds give different
  ## pi. Algorithm 1: "Generate a random one-to-one mapping pi : I -> [n]".
  pi_vec <- sample.int(n)           # pi_vec[i] = pi(i)
  pi_inv <- integer(n)              # pi_inv[pi(i)] = i  (the inverse map)
  pi_inv[pi_vec] <- seq_len(n)

  ## Split 1..n into consecutive blocks of size B; any tail of < B leftover
  ## indices is held fixed (it cannot be cyclically shifted within a full
  ## block).
  ## Algorithm 1, case 2 is i <= (K+1) floor(n/(K+1)): these are the indices
  ## that are shifted. Case 1, i > (K+1) floor(n/(K+1)), is the tail that
  ## psi~_k(i) = i leaves alone.
  nb <- n %/% B                     # number of full blocks
  in_block <- seq_len(nb * B)       # the indices that live in a full block
  pos0 <- (in_block - 1L) %% B      # 0-based position within the block
  blk_start <- in_block - pos0      # 1-based index of the block's first element

  ## Each non-identity element k cyclically shifts every block by k positions;
  ## element 0 is the identity. Building all powers of one generator guarantees
  ## the set is a cyclic group of order B (closed under composition).
  perms <- vector("list", B)
  for (k in 0:K) {
    psit <- seq_len(n)              # psi-tilde_k image vector; identity
    if (k > 0L && nb > 0L) {
      ## Algorithm 1's display: psi~_k(i) = i + k if the residue is <= K+1-k,
      ## else i - (K+1-k). With the 1-based residue p = pos0 + 1 that is the
      ## cyclic shift p -> ((p - 1 + k) mod (K+1)) + 1 of every block.
      new_pos0 <- (pos0 + k) %% B   # shifted position within the block
      psit[in_block] <- blk_start + new_pos0
    }
    ## Conjugate the block shift by the relabelling: psi_k(i) = pi^{-1}(
    ## psi-tilde_k( pi(i) ) )
    perms[[k + 1L]] <- pi_inv[psit[pi_vec]]
  }
  attr(perms, "block_size") <- B
  perms
}

## ---- internal RNG-state helpers ------------------------------------------
## These let a seeded routine draw random numbers reproducibly without
## clobbering the caller's random stream: snapshot .Random.seed, reseed, do the
## work, then restore the snapshot on exit.

#' Snapshot the current global RNG state (or NULL if none has been
#' initialised).
#' @details RNG hygiene around the seeded draws of Algorithm 1 (and of every
#'   other seeded step); not itself a paper step.
#' @keywords internal
#' @noRd
.save_seed <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
}

#' Restore a global RNG state previously captured by `.save_seed()`. A NULL
#' snapshot means the RNG was uninitialised, so we remove the seed again.
#' @details RNG hygiene around the seeded draws of Algorithm 1; not itself a
#'   paper step.
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
