## Internal computational core for the invariant permutation test.
## Not exported.

#' Residual maker via QR (rank robust) -- REFERENCE IMPLEMENTATION
#'
#' Returns the residuals of `V` after projecting onto the column space of `M`,
#' i.e. (I - P_M) V where P_M is the orthogonal projector onto col(M). Uses a
#' pivoted QR so that rank-deficient `M` (e.g. a duplicated intercept after
#' stacking X with a permuted copy of X) is handled correctly.
#'
#' **Not on the fit path.** Nothing in R/ calls this. `.ipt_prepare()`
#' computes the same residuals through an eigendecomposition pseudo-inverse of
#' the 2p x 2p Gram matrix: that never materializes the N x 2p stack, and it
#' replaces a QR of that stack (O(N (2p)^2)) with one O(N p^2) cross product
#' plus an O(p^3) eigendecomposition that does not grow with N. This function
#' survives as the obvious, obviously correct formulation that
#' `tests/lower-level-tests/test-projection.R` pins against an independent SVD
#' projector; keep the two in agreement, and do not delete this without
#' relocating that test.
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
#' profiling. This loop returns the same numbers (max is a comparison
#' reduction: no arithmetic, so no reassociation to worry about) without
#' copying the matrix. Names are dropped -- every caller indexes by position.
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
#' Maps a matrix of integer coordinates (each column taking values in 1..max)
#' to a single numeric mixed-radix code, so cells can be matched with
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

#' Within-cell slot index: the position of each observation inside its cell
#'
#' The cell-code machinery keys an observation by its cluster coordinates, and
#' those coordinates identify an observation uniquely only when each cell
#' holds one observation. Designs with replication inside a cell (two-way
#' layouts) therefore need a third coordinate: the slot l = 1..ell_ij that an
#' observation occupies inside its cell. `mwperm_layout()` permutes within a
#' cell, so its slot is a within-cell RANK; `mwperm_irregular()` holds the
#' slot fixed across cells, so its slot must be the `rep` LEVEL itself (the
#' same period in every cell -- see `.irregular_design()`), and it calls this
#' function only for the `rep = NULL` case, where the rank is the level.
#'
#' The ordering is the one `mwperm_layout()` has always used: by `rep` where
#' supplied (as a factor, so labels of any type order consistently), by order
#' of appearance otherwise, with ties broken by position. Equivalent to
#' `rank(ties.method = "first")` within each cell, computed as one stable
#' sort.
#'
#' @param cell integer vector of dense 1-based cell ids, one per observation.
#' @param rep optional within-cell replication identifier; `NULL` means use
#'   the order of appearance.
#' @param ncell number of cells (the maximum cell id).
#' @return integer vector of slot indices, running 1..ell inside each cell.
#' @keywords internal
#' @noRd
.within_cell_slot <- function(cell, rep = NULL, ncell = max(cell)) {
  N <- length(cell)
  slot <- integer(N)
  ord_key <- if (is.null(rep)) seq_len(N) else as.numeric(factor(rep))
  o <- order(cell, ord_key)
  slot[o] <- sequence(tabulate(cell, nbins = ncell))
  slot
}

#' Apply one group element to the rows of a data matrix (or vector)
#'
#' The engine's group elements are *signed gathers*: a two-slot list
#' `list(g = <integer gather vector or NULL>, s = <+/-1 vector or NULL>)`,
#' applied as `M[g, ] * s`. `NULL` in a slot means the identity for that slot,
#' so a permutation element carries `s = NULL`, a sign-flip element carries
#' `g = NULL`, and a combined permute-and-flip element carries both. A bare
#' integer vector is accepted as a pure gather and is applied by exactly the
#' indexing expression the permutation front ends have always used, so those
#' front ends need no change and their seeded output is bit-identical.
#'
#' The representation is two independent slots rather than a type tag with
#' one payload for a specific reason: the algebra that `.ipt_prepare()` relies
#' on -- X_k = Q_k X for an ORTHOGONAL Q_k, so that X_k' X_k = X' X -- holds
#' for a permutation matrix, for a diagonal +/-1 matrix, and for their product
#' alike. Any element of the group generated by both kinds is therefore a
#' single `(g, s)` pair, and a future design that permutes clusters *and*
#' flips signs (the two invariances are complementary, not nested) composes
#' its elements without touching the engine or this function. A tagged union
#' would have forced that design to add a third tag and a third branch here.
#'
#' Multiplication by +/-1 is exact in floating point, so a sign-flip element
#' introduces no rounding; the sign is applied AFTER the gather so that
#' `.apply_op(op, M)` reads as "gather the rows, then flip them" and the
#' composition rule is `(g2, s2) o (g1, s1) = (g1[g2], s1[g2] * s2)`.
#'
#' @param op a group element: an integer gather vector, or a list with slots
#'   `g` (integer gather vector or `NULL`) and `s` (numeric +/-1 vector of
#'   length N or `NULL`).
#' @param M a numeric matrix with N rows, or a numeric vector of length N.
#' @return `M` with the element applied to its rows, same shape as `M`.
#' @keywords internal
#' @noRd
.apply_op <- function(op, M) {
  if (!is.list(op))                    # bare gather vector: the historical form
    return(if (is.matrix(M)) M[op, , drop = FALSE] else M[op])
  out <- M
  if (!is.null(op$g))
    out <- if (is.matrix(out)) out[op$g, , drop = FALSE] else out[op$g]
  if (!is.null(op$s)) out <- out * op$s   # recycles down the rows (column-major)
  out
}

#' Precompute the beta-independent pieces of Procedure 1
#'
#' Implements the expensive, beta-independent part of Procedure 1 of Guo,
#' Toulis and Wang (2026). For each non-identity group element k = 1..K it
#' partials out the nuisance design by residualizing on M_k = [X | X_k] (the
#' Frisch-Waugh-Lovell step, equivalent to projecting with the orthonormal V_k
#' that satisfies V_k' X = V_k' X_k = 0) and stores the small d-dimensional
#' cross products needed to evaluate the test statistic at *any* null value
#' beta = b. Because residualization is linear in the outcome, the statistics
#' a_k(b) = || D' V_k V_k' (y - D b) ||, b_k(b) = || D' V_k V_k' (y - D b)_k
#' || are affine in `b` once the residualized covariate Dr = V_k V_k' D and
#' its inner products with `y`, `y_k`, `D` and `D_k` are known. Caching them
#' here lets `.invert_ci()` sweep many candidate values of `b` without redoing
#' a single QR decomposition.
#'
#' **The group-element contract (widened for the sign-flip design).** An
#' element is anything `.apply_op()` accepts: a bare integer gather vector
#' (every permutation front end) or a signed gather `list(g = , s = )`
#' (`mwperm_dyadic_het()`, and any future permute-and-flip design). Here
#' `X_k`, `y_k` and `D_k` are `.apply_op(op, X)`, `.apply_op(op, y)` and
#' `.apply_op(op, D)`; nothing else in this function depends on which kind of
#' element it is given, because everything it exploits -- X_k' X_k = X' X, and
#' the invariance identity Dr' y_k = Dr' y when op leaves Dr fixed -- holds
#' for any orthogonal row action, and a signed permutation is orthogonal. The
#' affine-in-`b` structure the cache encodes is likewise unchanged: for the
#' sign-flip element S_k, b_k(b) = |Dr' S_k y - Dr' S_k D b| with Dr the
#' residual of D on [X | S_k X]. So `.ipt_eval()`, `.pval_matrix()`, the
#' exact confidence set and `.invert_ci()` all run on a sign-flip prep object
#' without knowing it is one.
#'
#' @param y numeric outcome, length N (the *unshifted* outcome).
#' @param D numeric N x d matrix of covariate(s) of interest.
#' @param X numeric N x p nuisance design (intercept already included).
#' @param obs_perms list of length K+1 of group elements over the N
#'   observations, each either an integer gather-vector or a signed gather
#'   `list(g, s)` (see `.apply_op()`); element 1 must be the identity (the
#'   engine never applies it, but the K+1 count is read from the list).
#' @param need_perm_D logical; if `FALSE` the permuted-`D` cross product `W`
#'   (needed only for confidence-interval inversion and non-zero nulls) is not
#'   formed. It is one inner product per permutation, so the saving is small;
#'   the field is left `NULL`-valued so a prep object built without it cannot
#'   silently be used for inversion.
#' @param n_cores number of workers for the per-permutation QR loop (the K
#'   independent factorizations); 1 = serial. The loop uses no RNG and each
#'   iteration writes an independent slice, so any schedule gives identical
#'   output.
#' @param degenerate what to do when the residualized `D` is numerically zero
#'   for some permutation k, i.e. no identifying variation in `d` survives the
#'   projection onto the orthogonal complement of `[X | X_k]`. `"stop"` (the
#'   default) raises an error naming the design: that slice's statistic is
#'   pure rounding noise, so the comparison a_k vs b_k is decided by float
#'   error rather than by the data, and silently returning a number would be
#'   worse than failing. `"zero"` restores the exact-arithmetic answer instead
#'   (a_k = b_k = 0, hence p = 1 through the minorization) and is used by the
#'   engine when it has ALREADY established that beta is unidentified -- `d`
#'   constant or collinear with `x` -- and warned about it. In that case p = 1
#'   is the correct answer, not a failure. (Procedure 1 defines the
#'   per-permutation case the same way -- a_k = b_k = 0, hence p = 1 -- so
#'   stopping is a deliberate refusal to report a value that floating point
#'   cannot certify, not a gap in the arithmetic; both branches and the
#'   relative threshold are pinned by `tests/lower-level-tests/test-pvalue.R`.)
#' @param design short label for the calling design, used in that error.
#' @return a list (the "prep" object) consumed by `.ipt_eval()`.
#' @keywords internal
#' @noRd
.ipt_prepare <- function(y, D, X, obs_perms, need_perm_D = TRUE, n_cores = 1L,
                         cl = NULL, degenerate = c("stop", "zero"),
                         design = "this design") {
  degenerate <- match.arg(degenerate)
  y <- as.numeric(y)
  D <- as.matrix(D)
  X <- as.matrix(X)
  Kp1 <- length(obs_perms)             # group order (incl. identity)
  K <- Kp1 - 1L                        # number of non-identity permutations
  if (K < 1L) stop("Need at least one non-identity permutation.", call. = FALSE)
  d <- ncol(D)                         # number of coefficients of interest
  ## Degeneracy floor: when D lies numerically inside span[X | X_k], the
  ## residualized cross products are pure rounding noise, and in exact
  ## arithmetic that slice has a_k = b_k = 0. Relative tolerance in the qr()
  ## league (1e-8 on the Frobenius norm). Two situations reach it and they
  ## deserve different answers, which is what `degenerate` selects. If beta is
  ## unidentified from the start (d constant, or collinear with x) then EVERY
  ## slice is degenerate, p = 1 is the exact answer, and zeroing the slice is
  ## right -- the engine detects that case up front, warns, and asks for
  ## "zero". Otherwise beta IS identified in the data and a degenerate slice
  ## means this particular permutation annihilated d: the statistic for that
  ## permutation would be decided by float noise, so the fit stops instead.
  degen_tol2 <- 1e-16 * sum(D * D)

  ## Every group element acts on the rows of X by an ORTHOGONAL matrix Q_k.
  ## For the permutation designs Q_k is a permutation matrix: every gather
  ## vector is a bijection (the complete-array builders match a permuted cell
  ## code back to a unique row; Procedure 2 permutes within bicliques that are
  ## disjoint in both margins). For the sign-flip design Q_k is a diagonal
  ## +/-1 matrix, and for a combined element it is their product; all three
  ## satisfy Q_k'Q_k = I. So X_k = Q_k X, and the Gram matrix of the stacked
  ## design M_k = [X | X_k] is
  ##     M_k'M_k = [ X'X    X'X_k ]
  ##               [ X_k'X   X'X  ]
  ## -- the lower-right block is X'Q_k'Q_k X = X'X, the SAME matrix as the
  ## upper-left, for every k. Only the cross block depends on the element.
  ## That is the whole optimization: X'X is formed once below, and the per-
  ## element cost drops to one crossprod (level-3 BLAS) in place of a QR
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

  ## One slice of cached cross products per non-identity element k.
  ## Naming: Dr = D residualized on M_k = [X | X_k]; a trailing `g` is the
  ## transformed (permuted and/or sign-flipped) copy of a vector, obtained
  ## through .apply_op() so that the same code serves both group kinds.
  one_k <- function(k) {
    op <- obs_perms[[k + 1L]]          # the k-th non-identity element
    ## Only D is residualized. V_k V_k' = I - P_{M_k} is symmetric and
    ## idempotent, so D' V_k V_k' z = (V_k V_k' D)' z = Dr' z for any z: the
    ## permuted and unpermuted outcomes enter as plain inner products against
    ## Dr. This is the statistic as Procedure 1 writes it,
    ## a_k(b) = ||Dr'(y - D b)||.
    Dr <- if (p == 0L) D else {
      Xg <- .apply_op(op, X)                          # X_k = Q_k X
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
    if (sum(Dr * Dr) <= degen_tol2) {  # degenerate slice
      if (degenerate == "stop")
        stop(sprintf(paste0("No identifying variation in `d` survives ",
                            "permutation %d of the %s design: the ",
                            "residualized `d` is numerically zero after ",
                            "projecting out the nuisance design and its ",
                            "permuted copy, so that permutation's test ",
                            "statistic would be rounding noise rather than ",
                            "data. This usually means `d` is (nearly) a ",
                            "linear combination of `x` and its permuted ",
                            "copy. Drop the redundant nuisance covariates, ",
                            "or test a `d` with independent variation."),
                     k, design), call. = FALSE)
      return(list(u = matrix(0, d, 1L), v = matrix(0, d, 1L),
                  M = matrix(0, d, d),
                  W = if (need_perm_D) matrix(0, d, d)))
    }
    ## Invariant slice. When the element leaves the residualized regressor
    ## bitwise unchanged -- Q_k Dr == Dr -- the transformed statistic EQUALS
    ## the untransformed one in exact arithmetic: Q_k is orthogonal, so
    ##     Dr' Q_k y = (Q_k' Dr)' y = (Q_k^{-1} Dr)' y = Dr' y,
    ## (for a permutation, sum_i Dr_i y_g(i) = sum_m Dr_g^-1(m) y_m) and
    ## likewise Dr' D_k = Dr' Dr. Recomputing them instead sums the same
    ## terms in a different order, so v and W come back differing from u and M
    ## by ~1e-13 of pure rounding noise, and that noise -- not the data --
    ## then decides the a_j <= b_k comparisons. This is the no-power case
    ## Section 6.3 warns about (`d` constant within cell, permuted within
    ## cell): the exact answer is a_k == b_k for every k, hence p = 1. Assert
    ## the identity rather than recompute it, so the answer is the same on
    ## every platform's BLAS. Non-degenerate slices never take this branch,
    ## and a non-identity sign flip cannot (it negates some row of Dr).
    uu <- crossprod(Dr, y)                         # = Dr' y
    MM <- crossprod(Dr)                            # = Dr' Dr
    if (identical(.apply_op(op, Dr), Dr))
      return(list(u = uu, v = uu, M = MM,
                  W = if (need_perm_D) MM))
    list(u = uu,
         v = crossprod(Dr, .apply_op(op, y)),      # = Dr' y_k
         M = MM,
         W = if (need_perm_D)          # only for CI / non-zero null
           crossprod(Dr, .apply_op(op, D)))
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
#' from `.ipt_prepare()`. Returns the minorized randomization p-value (1 +
#' sum_k 1{ min_j a_j(b) <= b_k(b) }) / (K + 1).
#'
#' @param prep a prep object from `.ipt_prepare()`.
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
#' Thin wrapper combining `.ipt_prepare()` and `.ipt_eval()` to evaluate the
#' test at beta = 0 for an already-shifted outcome `y`. Kept for direct use
#' and testing; the engine uses prepare/eval separately so the QR work is
#' shared across the confidence-interval search.
#'
#' @inheritParams .ipt_prepare
#' @return list with `pvalue`, and diagnostic vectors `a`, `b`.
#' @keywords internal
#' @noRd
.ipt_pvalue <- function(y, D, X, obs_perms) {
  prep <- .ipt_prepare(y, D, X, obs_perms, need_perm_D = FALSE)
  .ipt_eval(prep, rep(0, prep$d))
}

#' Is a position table over the whole cell index space worth allocating?
#'
#' The two gather-vector builders translate a permuted cell code back to an
#' observation index either through an integer table indexed by the code --
#' one O(N) gather per group element -- or through `match()`. The table is
#' faster but its size is the whole mixed-radix index space, which for sparse
#' ids can dwarf the data: 300 x 300 x 700 levels is 63 million cells. It
#' used to be allocated up to 2^26 entries (268 MB) regardless of N. Two caps
#' now apply: at most 2^24 entries (64 MB), and at most 64 entries per
#' observation, so the table is never more than a small multiple of the data
#' it indexes. The output is identical on either branch (asserted by
#' `tests/lower-level-tests/test-obsperms.R`, which forces both).
#'
#' @param n_cells size of the index space (product of the radices).
#' @param N number of observations.
#' @keywords internal
#' @noRd
.use_pos_table <- function(n_cells, N)
  n_cells <= 2^24 && n_cells <= 64 * as.double(N)

#' Build observation-level permutations from per-dimension permutation groups
#'
#' Given the integer cluster id of every observation along each clustering
#' dimension, and a list of per-dimension permutation groups (each a list of
#' K+1 image vectors), produce the K+1 observation-level gather vectors. A
#' coordinate dimension can be held fixed by passing `NULL` for its group.
#'
#' Every returned gather vector is checked to be a genuine permutation of
#' `seq_len(N)` before it is handed back (see `.assert_bijection()`): the
#' cells are keyed by a mixed-radix code and translated back to observation
#' indices, and a code shared by two observations would make that translation
#' many-to-one, silently computing the statistic on duplicated rows. Duplicate
#' cell codes are therefore rejected up front as well.
#'
#' @param coords integer matrix N x C of cluster ids (1-based, dense).
#' @param groups list of length C; each element is either `NULL` (dimension
#'   held fixed) or a list of K+1 image vectors permuting that dimension's
#'   ids.
#' @param design short label for the calling design, used in error messages.
#' @param front_end the front end a user should reach for instead, named in
#'   the duplicate-cell error.
#' @param pos_table `NULL` (the default) lets `.use_pos_table()` decide from
#'   the size of the index space whether permuted cell codes are translated
#'   back through a position table or through `match()`; `TRUE`/`FALSE`
#'   forces one branch. Both give identical gather vectors -- the tests force
#'   each and compare -- so this exists only to make that assertion possible.
#' @return list of length K+1 of integer gather-vectors over observations.
#' @keywords internal
#' @noRd
.build_obs_perms <- function(coords, groups, design = "this design",
                             front_end = paste("mwperm_layout() or",
                                               "mwperm_missing()"),
                             pos_table = NULL) {
  coords <- as.matrix(coords)
  C <- ncol(coords)                    # number of clustering dimensions
  ## Mixed-radix bases, computed ONCE for all K+1 elements. Each image vector
  ## is a bijection of its dimension's ids, so a permutation cannot change a
  ## column maximum: the permuted coordinates have the same radix as `coords`.
  ## The old code recomputed it inside .cell_code on every element (and twice
  ## more for the `pos` table below) -- see .col_max for why that was costly.
  radix <- .col_max(coords)
  orig_code <- .cell_code(coords, radix = radix)   # cell code of each obs
  ## Cells must be unique. Below, a permuted cell code is translated back to an
  ## observation index by match() / a position table, both of which return the
  ## FIRST observation carrying that code. If two observations shared a cell,
  ## the translation would be many-to-one: the "gather vector" would repeat one
  ## row and drop another, and the statistic would be computed on duplicated
  ## data with no error raised anywhere. Reject that here rather than
  ## discovering it as a bijection failure per element below.
  dup <- anyDuplicated(orig_code)
  if (dup > 0L)
    stop(sprintf(paste0("%s requires exactly one observation per cell, but ",
                        "cell (%s) appears more than once (observation %d). ",
                        "Repeated cells are within-cell replication or ",
                        "repeated time periods -- use %s."),
                 design,
                 paste(coords[dup, ], collapse = ", "), dup, front_end),
         call. = FALSE)
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
  if (if (is.null(pos_table)) .use_pos_table(n_cells, nrow(coords))
      else isTRUE(pos_table)) {
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
  .assert_bijection(obs_perms, nrow(coords), design)
  obs_perms
}

#' Assert that every gather vector is a permutation of seq_len(N)
#'
#' The whole method rests on the permuted data being a relabelling of the
#' observed data: Procedure 1 needs X_k = Pi_k X for a permutation matrix
#' Pi_k, and .ipt_prepare() exploits that identity directly (it reuses X'X as
#' the lower-right Gram block, which is only correct for a bijection). A
#' gather vector that repeated an index would give a wrong statistic silently
#' -- no NA, no warning, just numbers computed on duplicated rows. This is
#' cheap next to a single permutation's linear algebra (O(N) per element
#' against O(N p^2)), so it runs unconditionally.
#'
#' @param obs_perms list of integer gather-vectors.
#' @param N expected length (the number of observations).
#' @param design short label for the calling design, used in the message.
#' @keywords internal
#' @noRd
.assert_bijection <- function(obs_perms, N, design = "this design") {
  for (k in seq_along(obs_perms)) {
    g <- obs_perms[[k]]
    if (length(g) != N || anyNA(g) || anyDuplicated(g))
      stop(sprintf(paste0("Internal error: permutation %d of the %s ",
                          "permutation group is not a bijection of the %d ",
                          "observations (length %d, %d distinct). The ",
                          "permuted data would then repeat some rows and ",
                          "drop others, so the test statistic would be ",
                          "wrong. Please report this with a reproducible ",
                          "example."),
                   k, design, N, length(g), length(unique(g))),
           call. = FALSE)
  }
  invisible(NULL)
}
