#' Construct a random block-cyclic permutation group (Algorithm 1)
#'
#' Builds a set of \code{K + 1} permutations of \code{seq_len(n)} that form a
#' cyclic algebraic group, following Algorithm 1 of Guo, Toulis and Wang (2026)
#' (adapted from Wen, Wang and Wang, 2025). A random one-to-one relabelling
#' \eqn{\pi} is drawn, the index set is split into consecutive blocks of size
#' \code{K + 1}, and each non-identity element cyclically shifts every block by
#' \eqn{k} positions. Because the elements are powers of a single generator,
#' the returned set is closed under composition (a cyclic group of order
#' \code{K + 1}).
#'
#' The randomness in the relabelling is what makes the resulting test a
#' \emph{random} invariant test: different seeds yield slightly different
#' p-values. Aggregating several runs (see the \code{n_reps} argument of the
#' \code{mwperm_*} functions) by taking the median p-value is recommended.
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
#' Guo, F. R., Toulis, P. and Wang, Y. (2026). Permutation inference under
#' multi-way clustering and missing data.
#'
#' Wen, K., Wang, T. and Wang, Y. (2025). Residual permutation test for
#' regression coefficient testing. \emph{The Annals of Statistics} 53(2),
#' 724--748.
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
  B <- K + 1L                      # block size
  if (B > n) {
    stop(sprintf(
      "K + 1 = %d exceeds n = %d: no non-trivial permutation exists. Use K <= n - 1.",
      B, n), call. = FALSE)
  }

  if (!is.null(seed)) {
    old <- .save_seed()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }

  ## random relabelling pi and its inverse
  pi_vec <- sample.int(n)           # pi_vec[i] = pi(i)
  pi_inv <- integer(n)
  pi_inv[pi_vec] <- seq_len(n)

  nb <- n %/% B                     # number of full blocks
  in_block <- seq_len(nb * B)       # indices that live in a full block
  ## position within block (0-based) and block start (1-based) for in-block idx
  pos0 <- (in_block - 1L) %% B
  blk_start <- in_block - pos0      # first index of each element's block

  perms <- vector("list", B)
  for (k in 0:K) {
    psit <- seq_len(n)              # psi-tilde_k as image vector; default identity
    if (k > 0L && nb > 0L) {
      new_pos0 <- (pos0 + k) %% B
      psit[in_block] <- blk_start + new_pos0
    }
    ## psi_k(i) = pi^{-1}( psi-tilde_k( pi(i) ) )
    perms[[k + 1L]] <- pi_inv[psit[pi_vec]]
  }
  attr(perms, "block_size") <- B
  perms
}

## ---- internal RNG-state helpers ------------------------------------------

.save_seed <- function() {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    NULL
  }
}

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
