## Internal computational core for the invariant permutation test.
## Not exported.

#' Residual maker via QR (rank robust) -- REFERENCE IMPLEMENTATION
#'
#' Returns the residuals of `V` after projecting onto the column space of `M`,
#' i.e. (I - P_M) V where P_M is the orthogonal projector onto col(M). Uses a
#' pivoted QR so that rank-deficient `M` (e.g. a duplicated intercept after
#' stacking X with a permuted copy of X) is handled correctly.
#'
#' **Not on the fit path.** Nothing in R/ calls this. `.ipt_prepare()` computes
#' the same residuals through an eigendecomposition pseudo-inverse of the
#' 2p x 2p Gram matrix: that never materializes the N x 2p stack, and it
#' replaces a QR of that stack (O(N (2p)^2)) with one O(N p^2) cross product
#' plus an O(p^3) eigendecomposition that does not grow with N. This function
#' survives as the obvious, obviously correct formulation that
#' `tests/test-projection.R` pins against an independent SVD projector; keep
#' the two in agreement, and do not delete this without relocating that test.
#'
#' @param M numeric matrix, N x q.
#' @param V numeric vector or N x d matrix.
#' @return residuals, same shape as `V`.
#' @keywords internal
#' @noRd
.residualize <- function(M, V) {
  qr.resid(qr(M), V)
}

#' Column maxima of a numeric matrix
#'
#' `apply(m, 2L, max)` routes the whole matrix through `aperm()` and a list
#' split, which showed up as ~a quarter of the permutation builder's time in
#' profiling. This loop returns the same numbers (max is a
#' comparison reduction: no arithmetic, so no reassociation to worry about)
#' without copying the matrix. Names are dropped -- every caller indexes by
#' position.
#'
#' @param m numeric matrix.
#' @return numeric vector of length `ncol(m)`.
#' @keywords internal
#' @noRd
.col_max <- function(m) {
  out <- numeric(ncol(m))
  for (c in seq_len(ncol(m))) out[c] <- max(m[, c])
  out
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
.cell_code <- function(coords, radix = NULL) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  ## Per-dimension base = number of clusters. `radix` may be supplied by a
  ## caller that already knows it (see .build_obs_perms, which reuses one
  ## radix for all K+1 permuted copies); the column maxima are the same
  ## numbers either way, so the codes are unchanged.
  if (is.null(radix)) radix <- .col_max(coords)
  ## Guard exactness: the code ranges over 0..prod(radix)-1, which must fit in a
  ## double without rounding (integers are exact only up to 2^53). For realistic
  ## cluster counts this is never an issue, but error loudly rather than risk a
  ## silent cell-code collision on a pathologically large design.
  if (prod(radix) > 2^53)
    stop(sprintf(paste0("Cluster index space (%.3g cells) exceeds 2^53; ",
                        "the mixed-radix cell encoding would not be exact. ",
                        "Too many clusters for this design."), prod(radix)),
         call. = FALSE)
  code <- coords[, 1L] - 1             # least-significant digit (0-based)
  if (ncol(coords) > 1L) {
    mult <- 1                          # place value of the current digit
    for (c in 2L:ncol(coords)) {
      mult <- mult * radix[c - 1L]     # advance to the next, higher base
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
#' are affine in `b` once the residualized covariate Dr = V_k V_k' D and its
#' inner products with `y`, `y_k`, `D` and `D_k` are known. Caching them here
#' lets \code{\link{.invert_ci}} sweep many candidate values of `b` without
#' redoing a single QR decomposition.
#'
#' @param y numeric outcome, length N (the *unshifted* outcome).
#' @param D numeric N x d matrix of covariate(s) of interest.
#' @param X numeric N x p nuisance design (intercept already included).
#' @param obs_perms list of length K+1 of integer gather-vectors over the N
#'   observations; element 1 must be the identity `seq_len(N)`.
#' @param need_perm_D logical; if `FALSE` the permuted-`D` cross product `W`
#'   (needed only for confidence-interval inversion and non-zero nulls) is not
#'   formed. It is one inner product per permutation, so the saving is small;
#'   the field is left `NULL`-valued so a prep object built without it cannot
#'   silently be used for inversion.
#' @param n_cores number of workers for the per-permutation QR loop (the K
#'   independent factorizations); 1 = serial. The loop uses no RNG and each
#'   iteration writes an independent slice, so any schedule gives identical
#'   output.
#' @return a list (the "prep" object) consumed by \code{\link{.ipt_eval}}.
#' @keywords internal
#' @noRd
.ipt_prepare <- function(y, D, X, obs_perms, need_perm_D = TRUE, n_cores = 1L,
                         cl = NULL) {
  y <- as.numeric(y)
  D <- as.matrix(D)
  X <- as.matrix(X)
  Kp1 <- length(obs_perms)             # group order (incl. identity)
  K <- Kp1 - 1L                        # number of non-identity permutations
  if (K < 1L) stop("Need at least one non-identity permutation.", call. = FALSE)
  d <- ncol(D)                         # number of coefficients of interest
  ## Degeneracy floor: when D lies numerically inside span[X | X_k], the
  ## residualized cross products are pure rounding noise. In exact arithmetic
  ## that slice has a_k = 0 -- forcing p = 1 through the minorization -- so
  ## the slice is zeroed to restore the exact answer rather than letting
  ## ~1e-16 noise decide the a/b comparisons. Relative tolerance
  ## in the qr() league (1e-8 on the Frobenius norm).
  degen_tol2 <- 1e-16 * sum(D * D)

  ## Every observation gather-vector is a bijection, in every design (the
  ## complete-array builders match a permuted cell code back to a unique row;
  ## Procedure 2 permutes within bicliques that are disjoint in both margins).
  ## So X_k = Pi_k X for a permutation matrix Pi_k, and the Gram matrix of the
  ## stacked design M_k = [X | X_k] is
  ##     M_k'M_k = [ X'X    X'X_k ]
  ##               [ X_k'X   X'X  ]
  ## -- the lower-right block is X'Pi_k'Pi_k X = X'X, the SAME matrix as the
  ## upper-left, for every k. Only the cross block depends on the permutation.
  ## That is the whole optimization: X'X is formed once below, and the per-
  ## permutation cost drops to one crossprod (level-3 BLAS) in place of a QR
  ## of the N x 2p stack.
  p <- ncol(X)                         # nuisance columns (incl. intercept)
  A <- crossprod(X)                    # shared by every permutation
  XtD <- crossprod(X, D)               # likewise
  ## X' held explicitly so the cross block below can be written as a plain
  ## product. `crossprod(X, Xg)` and `tX %*% Xg` are the same mathematical
  ## sums over the same reduction index; they differ only in how they enter
  ## BLAS, through dgemm 'T' and dgemm 'N' respectively. Under the reference
  ## BLAS the 'T' kernel accumulates each entry in a single scalar, so its
  ## FMA chain is serial (~1.4 GFLOP/s here), while the 'N' kernel's inner
  ## loop is a contiguous length-p axpy with p independent accumulators --
  ## same arithmetic, ~1.5x the speed, on the term that is two thirds of
  ## this function. An optimized BLAS blocks and reassociates the reduction
  ## and may do so differently for the two kernels, so bitwise agreement
  ## between the two spellings is a property of the BLAS, not of the algebra:
  ## verified on the reference BLAS at 67 of 67 permutations, max |diff|
  ## 0.000e+00.
  tX <- t(X)

  ## One slice of cached cross products per non-identity permutation k.
  ## Naming: Dr = D residualized on M_k = [X | X_k]; a trailing [g] is the
  ## permuted (gathered) copy of a vector.
  one_k <- function(k) {
    g <- obs_perms[[k + 1L]]           # gather-vector of the k-th permutation
    ## Only D is residualized. V_k V_k' = I - P_{M_k} is symmetric and
    ## idempotent, so D' V_k V_k' z = (V_k V_k' D)' z = Dr' z for any z: the
    ## permuted and unpermuted outcomes enter as plain inner products against
    ## Dr. This is the statistic as Procedure 1 writes it,
    ## a_k(b) = ||Dr'(y - D b)||.
    Dr <- if (p == 0L) D else {
      Xg <- X[g, , drop = FALSE]                      # X_k = Pi_k X
      B <- tX %*% Xg                                  # = X'X_k, the only
                                                      #   O(N p^2) term
      ## M_k is rank-deficient by construction -- the permuted intercept
      ## duplicates the original, and panel time dummies are permuted among
      ## themselves -- so the projector goes through a pseudo-inverse. The
      ## eigenvalues of the Gram matrix are the SQUARED singular values of
      ## M_k, so the cut below at 1e-14 is the same relative tolerance on
      ## singular values (1e-7) that .lm.fit applies, squared.
      e <- eigen(rbind(cbind(A, B), cbind(t(B), A)), symmetric = TRUE)
      pos <- e$values > 1e-14 * e$values[1L]
      V <- e$vectors[, pos, drop = FALSE]
      cf <- V %*% ((1 / e$values[pos]) *
                     crossprod(V, rbind(XtD, crossprod(Xg, D))))
      D - X %*% cf[seq_len(p), , drop = FALSE] -
          Xg %*% cf[p + seq_len(p), , drop = FALSE]   # (I - P_{M_k}) D
    }
    if (sum(Dr * Dr) <= degen_tol2)    # degenerate slice: exact-arithmetic zero
      return(list(u = matrix(0, d, 1L), v = matrix(0, d, 1L),
                  M = matrix(0, d, d),
                  W = if (need_perm_D) matrix(0, d, d)))
    list(u = crossprod(Dr, y),                     # = Dr' y
         v = crossprod(Dr, y[g]),                  # = Dr' y_k
         M = crossprod(Dr),                        # = Dr' Dr
         W = if (need_perm_D)          # only for CI / non-zero null
           crossprod(Dr, D[g, , drop = FALSE]))
  }
  slices <- .plapply(seq_len(K), one_k, n_cores = n_cores, cl = cl)

  ## Assemble in k order (order-independent: each slice is self-contained).
  u <- matrix(0, d,
              K)                 # u[, k]   = Dr' yr    (a-statistic intercept)
  v <- matrix(0, d,
              K)                 # v[, k]   = Dr' ypr   (b-statistic intercept)
  M <- array(0, dim = c(d, d,
                        K))      # M[ , ,k] = Dr' Dr    (a-statistic slope in b)
  W <- array(0, dim = c(d, d,
                        K))      # W[ , ,k] = Dr' Dpr   (b-statistic slope in b)
  for (k in seq_len(K)) {
    u[, k]   <- slices[[k]]$u
    v[, k]   <- slices[[k]]$v
    M[, , k] <- slices[[k]]$M
    if (need_perm_D) W[, , k] <- slices[[k]]$W
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
  beta <- rep(as.numeric(beta),
              length.out = d)   # recycle scalar null to length d
  ## a[k], b[k]: the identity- and k-th-permutation residual norms at this beta,
  ## reconstructed from the cached cross products (affine in beta, no QR).
  if (d == 1L) {                       # scalar fast path: norms reduce to abs()
    a <- abs(prep$u[1L, ] - prep$M[1L, 1L, ] * beta)
    b <- abs(prep$v[1L, ] - prep$W[1L, 1L, ] * beta)
  } else {
    K <- prep$K
    a <- numeric(K)
    b <- numeric(K)
    for (k in seq_len(K)) {
      a[k] <- sqrt(sum((prep$u[, k] - prep$M[, , k] %*% beta)^2))
      b[k] <- sqrt(sum((prep$v[, k] - prep$W[, , k] %*% beta)^2))
    }
  }
  amin <- min(a)                       # minorizing identity statistic
  ## Minorized randomization p-value: fraction of permutations whose statistic
  ## is at least the (minorized) observed one, with the usual +1 correction.
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
  C <- ncol(coords)                    # number of clustering dimensions
  ## Mixed-radix bases, computed ONCE for all K+1 elements. Each image vector
  ## is a bijection of its dimension's ids, so a permutation cannot change a
  ## column maximum: the permuted coordinates have the same radix as `coords`.
  ## The old code recomputed it inside .cell_code on every element (and twice
  ## more for the `pos` table below) -- see .col_max for why that was costly.
  radix <- .col_max(coords)
  orig_code <- .cell_code(coords, radix = radix)   # cell code of each obs
  ## Determine the group order K+1 from the first dimension that is permuted
  ## (all non-NULL groups share the same order).
  Kp1 <- NULL
  for (g in groups) {
    if (!is.null(g)) {
      Kp1 <- length(g)
      break
    }
  }
  if (is.null(Kp1)) stop("At least one dimension must be permuted.",
                         call. = FALSE)

  ## Position table over the mixed-radix index space: translating permuted cell
  ## codes through it replaces match()'s per-k double hashing with one O(N)
  ## integer gather (identical output; ~2x builder). Only when
  ## the index space is cheap to allocate; huge sparse spaces keep match().
  n_cells <- prod(radix)
  pos <- NULL
  if (n_cells <= 2^26) {
    pos <- integer(n_cells)            # 0 = unobserved cell
    pos[orig_code + 1] <- seq_len(nrow(coords))
  }

  ## Digits of the ORIGINAL coordinates, as doubles, extracted once. The
  ## per-element loop below reads these instead of copying `coords` and
  ## re-coercing it on every k.
  digit <- vector("list", C)
  for (c in seq_len(C)) digit[[c]] <- as.double(coords[, c])

  obs_perms <- vector("list", Kp1)
  for (k in seq_len(Kp1)) {
    ## Cell code of the coordinates under the k-th element, accumulated digit
    ## by digit. This is .cell_code() inlined over the permuted coordinates:
    ## same bases, same place values, same left-to-right accumulation order,
    ## so the codes are bit-for-bit what .cell_code(mapped) would return --
    ## it just never materializes `mapped`.
    mult <- 1                          # place value of the current digit
    mapped_code <- NULL
    for (c in seq_len(C)) {
      gk <- groups[[c]]
      ## NULL group => dimension held fixed (this is how panel freezes time)
      dc <- if (is.null(gk)) digit[[c]] else as.double(gk[[k]])[coords[, c]]
      if (c == 1L) {
        mapped_code <- dc - 1          # least-significant digit (0-based)
      } else {
        mult <- mult * radix[c - 1L]   # advance to the next, higher base
        mapped_code <- mapped_code + mult * (dc - 1)
      }
    }
    ## Translate permuted coordinates back to observation indices; an NA (or a
    ## zero table entry) means the permutation reached an unobserved cell.
    g <- if (is.null(pos)) match(mapped_code, orig_code) else {
      gi <- pos[mapped_code + 1]
      gi[gi == 0L] <- NA_integer_
      gi
    }
    if (anyNA(g)) {
      stop("Permutation maps to an unobserved cell; the design is not a ",
           "complete array. Use mwperm_missing() or mwperm_layout().",
           call. = FALSE)
    }
    obs_perms[[k]] <- g                # observation gather-vector for element k
  }
  obs_perms
}
