## Internal computational core for the invariant permutation test.
## Not exported.

#' Residual maker via QR (rank robust)
#'
#' Returns the residuals of `V` after projecting onto the column space of `M`,
#' i.e. (I - P_M) V where P_M is the orthogonal projector onto col(M). Uses a
#' pivoted QR so that rank-deficient `M` (e.g. a duplicated intercept after
#' stacking X with a permuted copy of X) is handled correctly.
#'
#' @param M numeric matrix, N x q.
#' @param V numeric vector or N x d matrix.
#' @return residuals, same shape as `V`.
#' @keywords internal
#' @noRd
.residualize <- function(M, V) {
  qr.resid(qr(M), V)
}

#' Encode integer cluster coordinates as a unique numeric code
#'
#' Maps a matrix of integer coordinates (each column taking values in
#' 1..max) to a single numeric mixed-radix code, so cells can be matched with
#' `match()`. Uses doubles to avoid 32-bit integer overflow.
#'
#' @param coords integer matrix, one row per observation, one column per
#'   clustering coordinate. Values must be positive integers.
#' @return numeric vector of codes, one per row.
#' @keywords internal
#' @noRd
.cell_code <- function(coords) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  radix <- apply(coords, 2L, max)
  ## Guard exactness: the code ranges over 0..prod(radix)-1, which must fit in a
  ## double without rounding (integers are exact only up to 2^53). For realistic
  ## cluster counts this is never an issue, but error loudly rather than risk a
  ## silent cell-code collision on a pathologically large design.
  if (prod(radix) > 2^53)
    stop(sprintf(paste0("Cluster index space (%.3g cells) exceeds 2^53; the ",
                        "mixed-radix cell encoding would not be exact. Too many ",
                        "clusters for this design."), prod(radix)), call. = FALSE)
  code <- coords[, 1L] - 1
  if (ncol(coords) > 1L) {
    mult <- 1
    for (c in 2L:ncol(coords)) {
      mult <- mult * radix[c - 1L]
      code <- code + mult * (coords[, c] - 1)
    }
  }
  code
}

#' Precompute the beta-independent pieces of Procedure 1
#'
#' Implements the expensive, beta-independent part of Procedure 1 of Guo, Toulis
#' and Wang (2026). For each non-identity permutation k = 1..K it partials out
#' the nuisance design by residualizing on M_k = [X | X_k] (the
#' Frisch-Waugh-Lovell step, equivalent to projecting with the orthonormal V_k
#' that satisfies V_k' X = V_k' X_k = 0) and stores the small d-dimensional
#' cross products needed to evaluate the test statistic at *any* null value
#' \eqn{\beta = b}. Because residualization is linear in the outcome, the
#' statistics
#'   a_k(b) = || D' V_k V_k' (y - D b) ||,
#'   b_k(b) = || D' V_k V_k' (y - D b)_k ||
#' are affine in `b` once the projections of `y`, `y_k`, `D` and `D_k` onto V_k
#' are known. Caching them here lets \code{\link{.invert_ci}} sweep many
#' candidate values of `b` without redoing a single QR decomposition.
#'
#' @param y numeric outcome, length N (the *unshifted* outcome).
#' @param D numeric N x d matrix of covariate(s) of interest.
#' @param X numeric N x p nuisance design (intercept already included).
#' @param obs_perms list of length K+1 of integer gather-vectors over the N
#'   observations; element 1 must be the identity `seq_len(N)`.
#' @param need_perm_D logical; if `FALSE` the permuted-`D` projection (needed
#'   only for confidence-interval inversion) is skipped to save work.
#' @return a list (the "prep" object) consumed by \code{\link{.ipt_eval}}.
#' @keywords internal
#' @noRd
.ipt_prepare <- function(y, D, X, obs_perms, need_perm_D = TRUE) {
  y <- as.numeric(y)
  D <- as.matrix(D)
  X <- as.matrix(X)
  Kp1 <- length(obs_perms)
  K <- Kp1 - 1L
  if (K < 1L) stop("Need at least one non-identity permutation.", call. = FALSE)
  d <- ncol(D)

  u <- matrix(0, d, K)                 # u[, k] = Dr' yr      (a-intercept)
  v <- matrix(0, d, K)                 # v[, k] = Dr' ypr     (b-intercept)
  M <- array(0, dim = c(d, d, K))      # M[ , , k] = Dr' Dr   (a-slope)
  W <- array(0, dim = c(d, d, K))      # W[ , , k] = Dr' Dpr  (b-slope)
  for (k in seq_len(K)) {
    g   <- obs_perms[[k + 1L]]
    qMk <- qr(cbind(X, X[g, , drop = FALSE]))
    Dr  <- qr.resid(qMk, D)            # N x d
    yr  <- qr.resid(qMk, y)            # N
    ypr <- qr.resid(qMk, y[g])         # N (permuted outcome residualized)
    u[, k]   <- crossprod(Dr, yr)
    v[, k]   <- crossprod(Dr, ypr)
    M[, , k] <- crossprod(Dr)
    if (need_perm_D)
      W[, , k] <- crossprod(Dr, qr.resid(qMk, D[g, , drop = FALSE]))
  }
  list(u = u, v = v, M = M, W = W, K = K, Kp1 = Kp1, d = d,
       has_perm_D = need_perm_D)
}

#' Evaluate the Procedure 1 p-value at a null value beta = b
#'
#' Cheap: O(K d^2), no matrix factorizations. Uses the cached cross products
#' from \code{\link{.ipt_prepare}}. Returns the minorized randomization p-value
#'   (1 + sum_k 1{ min_j a_j(b) <= b_k(b) }) / (K + 1).
#'
#' @param prep a prep object from \code{\link{.ipt_prepare}}.
#' @param beta numeric null value(s); recycled to length `prep$d`.
#' @return list with `pvalue`, and diagnostic vectors `a`, `b`.
#' @keywords internal
#' @noRd
.ipt_eval <- function(prep, beta) {
  d <- prep$d
  beta <- rep(as.numeric(beta), length.out = d)
  if (d == 1L) {
    a <- abs(prep$u[1L, ] - prep$M[1L, 1L, ] * beta)
    b <- abs(prep$v[1L, ] - prep$W[1L, 1L, ] * beta)
  } else {
    K <- prep$K
    a <- numeric(K); b <- numeric(K)
    for (k in seq_len(K)) {
      a[k] <- sqrt(sum((prep$u[, k] - prep$M[, , k] %*% beta)^2))
      b[k] <- sqrt(sum((prep$v[, k] - prep$W[, , k] %*% beta)^2))
    }
  }
  amin <- min(a)
  list(pvalue = (1 + sum(b >= amin)) / prep$Kp1, a = a, b = b)
}

#' Core invariant permutation p-value (Procedure 1), convenience wrapper
#'
#' Thin wrapper combining \code{\link{.ipt_prepare}} and \code{\link{.ipt_eval}}
#' to evaluate the test at \eqn{\beta = 0} for an already-shifted outcome `y`.
#' Kept for direct use and testing; the engine uses prepare/eval separately so
#' the QR work is shared across the confidence-interval search.
#'
#' @inheritParams .ipt_prepare
#' @return list with `pvalue`, and diagnostic vectors `a`, `b`.
#' @keywords internal
#' @noRd
.ipt_pvalue <- function(y, D, X, obs_perms) {
  prep <- .ipt_prepare(y, D, X, obs_perms, need_perm_D = FALSE)
  .ipt_eval(prep, rep(0, prep$d))
}

#' Build observation-level permutations from per-dimension permutation groups
#'
#' Given the integer cluster id of every observation along each clustering
#' dimension, and a list of per-dimension permutation groups (each a list of
#' K+1 image vectors), produce the K+1 observation-level gather vectors. A
#' coordinate dimension can be held fixed by passing `NULL` for its group.
#'
#' @param coords integer matrix N x C of cluster ids (1-based, dense).
#' @param groups list of length C; each element is either `NULL` (dimension
#'   held fixed) or a list of K+1 image vectors permuting that dimension's ids.
#' @return list of length K+1 of integer gather-vectors over observations.
#' @keywords internal
#' @noRd
.build_obs_perms <- function(coords, groups) {
  coords <- as.matrix(coords)
  C <- ncol(coords)
  orig_code <- .cell_code(coords)
  ## determine K from the first non-NULL group
  Kp1 <- NULL
  for (g in groups) if (!is.null(g)) { Kp1 <- length(g); break }
  if (is.null(Kp1)) stop("At least one dimension must be permuted.", call. = FALSE)

  obs_perms <- vector("list", Kp1)
  for (k in seq_len(Kp1)) {
    mapped <- coords
    for (c in seq_len(C)) {
      gk <- groups[[c]]
      if (!is.null(gk)) {
        img <- gk[[k]]                 # image vector for this dim at element k
        mapped[, c] <- img[coords[, c]]
      }
    }
    mapped_code <- .cell_code(mapped)
    g <- match(mapped_code, orig_code)
    if (anyNA(g)) {
      stop("Permutation maps to an unobserved cell; the design is not a ",
           "complete array. Use mwperm_missing() or mwperm_layout().",
           call. = FALSE)
    }
    obs_perms[[k]] <- g
  }
  obs_perms
}
