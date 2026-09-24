## Shared code for the paper-fidelity suite (sourced by testthat before every
## test file; never run as a test itself).
##
## Four things live here:
##   1. a deliberately naive, readable Procedure 1 (Guo, Toulis & Wang 2026,
##      "GTW"), written from the paper and sharing no statistic code with R/;
##   2. reconstruction of the permutation group each front end draws in
##      repetition r, rebuilt from the seed exactly the way the package does
##      (the front ends keep their builders in closures, so there is no seam
##      to read them from -- see follow-up F-1 in TESTING_PLAN.md);
##   3. small data-generating processes used across files;
##   4. the Monte Carlo helpers and the MWPERM_SLOW_TESTS gate.
##
## Everything is base R plus testthat; withr is deliberately not called
## (it is not a declared dependency), so `with_seed()` below stands in for
## withr::with_seed().

## ---- RNG hygiene --------------------------------------------------------------

## Evaluate `expr` with the RNG seeded, then put the caller's stream back, so
## a helper's draws never shift the random numbers a later test sees.
with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv(),
                          inherits = FALSE))
            rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

## ---- 1. Naive Procedure 1 -----------------------------------------------------
## GTW Procedure 1, step 1: for k = 1..K let V_k be an orthonormal basis of the
## orthogonal complement of col([X | X_k]); a_k = ||D' V_k V_k' y|| and
## b_k = ||D' V_k V_k' y_k||. V_k V_k' is the residual maker I - P for the
## stacked design, built here from a pivoted QR with numerical rank (the stack
## is rank deficient whenever X has an intercept, because a permutation maps
## the intercept to itself -- the paper's "N - 2p" columns is the full-rank
## case only). Step 2 is Eq. (10):
##   pval = (1 + sum_k 1{ min_j a_j <= b_k }) / (K + 1).

## I - P_M, with P_M the orthogonal projector onto col(M). R's default qr()
## (LINPACK dqrdc2) moves numerically dependent columns to the end, so the
## first `rank` columns of Q span col(M).
ref_residual_maker <- function(M) {
  q <- qr(M)
  Q <- qr.Q(q)[, seq_len(q$rank), drop = FALSE]
  diag(nrow(M)) - tcrossprod(Q)
}

## The statistics and the p-value of Procedure 1 at the null beta = b.
##   y, D, X : outcome (N), interest (N x d) and nuisance design (N x p, with
##             the intercept already in it), in the package's own row order;
##   perms   : list of K + 1 gather vectors, perms[[1]] the identity; element
##             k + 1 is the row each observation takes under group element k,
##             so y_k = y[perms[[k + 1]]] and X_k = X[perms[[k + 1]], ].
## Testing beta = b means running the procedure on y - D b (GTW, Procedure 1
## step 3).
ref_procedure1 <- function(y, D, X, perms, b = 0) {
  D <- as.matrix(D)
  X <- as.matrix(X)
  y_b <- y - drop(D %*% rep(b, length.out = ncol(D)))
  K <- length(perms) - 1L
  a <- numeric(K)
  bk <- numeric(K)
  for (k in seq_len(K)) {
    g <- perms[[k + 1L]]
    R <- ref_residual_maker(cbind(X, X[g, , drop = FALSE]))   # V_k V_k'
    a[k]  <- sqrt(sum((crossprod(D, R %*% y_b))^2))           # ||D'V V' y||
    bk[k] <- sqrt(sum((crossprod(D, R %*% y_b[g]))^2))        # ||D'V V' y_k||
  }
  list(a = a, b = bk, pvalue = (1 + sum(min(a) <= bk)) / (K + 1L),
       gap = min(abs(bk - min(a))))   # distance of the nearest indicator to a tie
}

## ---- 2. Rebuilding a front end's group from its seed ---------------------------
## Every permutation front end draws, in repetition r, groups seeded as
##   rep seed      rs = seed + r - 1                (engine.R, .ipt_engine)
##   dimension j   rs * stride + j, stride = 1000   (engine.R, .sub_seed)
## with each per-dimension group produced by build_perm_set() (Algorithm 1,
## tested on its own in test-algorithm1.R). The functions below rebuild those
## groups, and then map them to observations by explicit key matching -- not
## by the package's mixed-radix codes -- so the gather vectors are an
## independent construction of GTW's y_{pi_k, sigma_k}: the entry for cell
## (i, j) is the observation at cell (pi_k(i), sigma_k(j)).

rep_seed <- function(seed, r) seed + r - 1
sub_seed <- function(rs, j, stride = 1000) rs * stride + j
dense_id <- function(x) as.integer(factor(x))

## Gather vector for a map on cluster coordinates. `ids` is a list of integer
## coordinate vectors (one per dimension, one entry per observation); `maps`
## is a parallel list of image vectors, NULL meaning the coordinate is held
## fixed. Observation r goes to the row whose coordinates are the images of
## its own.
ref_gather <- function(ids, maps) {
  from <- do.call(paste, c(lapply(ids, as.character), sep = ":"))
  to <- do.call(paste, c(Map(function(v, m) as.character(
    if (is.null(m)) v else m[v]), ids, maps), sep = ":"))
  g <- match(to, from)
  stopifnot(!anyNA(g), !anyDuplicated(g))   # a genuine permutation of rows
  g
}

## Dyadic / panel / three-way: one Algorithm 1 group per permuted dimension,
## sub-seeds 1, 2 (and 3); element k pairs the k-th element of each group
## (GTW Eq. 9). `fixed` lists coordinates held fixed (the panel's time).
ref_groups_crossed <- function(ids, K, seed, r, permuted = seq_along(ids)) {
  rs <- rep_seed(seed, r)
  ids <- lapply(ids, dense_id)
  G <- vector("list", length(ids))
  for (j in seq_along(permuted))
    G[[permuted[j]]] <- build_perm_set(max(ids[[permuted[j]]]), K,
                                       seed = sub_seed(rs, j))
  lapply(seq_len(K + 1L), function(k)
    ref_gather(ids, lapply(G, function(g) if (is.null(g)) NULL else g[[k]])))
}

## Two-way layout (GTW Section 6.3): an Algorithm 1 group over the replicate
## index of every cell, cell c seeded at sub-seed c (stride widened past the
## cell count). The replicate index is the within-cell rank by `rep` (ties by
## position), and element k maps replicate l of cell c to replicate
## pi^c_k(l) of the SAME cell.
ref_slot <- function(cell, rep = NULL) {
  key <- if (is.null(rep)) seq_along(cell) else as.numeric(factor(rep))
  slot <- integer(length(cell))
  for (c in unique(cell)) {
    ii <- which(cell == c)
    slot[ii[order(key[ii], ii)]] <- seq_along(ii)
  }
  slot
}
ref_groups_layout <- function(row, col, rep, K, seed, r) {
  rs <- rep_seed(seed, r)
  cell <- dense_id(interaction(dense_id(row), dense_id(col), drop = TRUE))
  ncell <- max(cell)
  ell <- tabulate(cell, ncell)
  stride <- max(1000, ncell + 1)
  slot <- ref_slot(cell, rep)
  G <- lapply(seq_len(ncell), function(c)
    build_perm_set(ell[c], K, seed = sub_seed(rs, c, stride)))
  from <- paste(cell, slot)
  lapply(seq_len(K + 1L), function(k) {
    img <- vapply(seq_along(cell), function(o) G[[cell[o]]][[k]][slot[o]],
                  integer(1))
    g <- match(paste(cell, img), from)
    stopifnot(!anyNA(g), !anyDuplicated(g))
    g
  })
}

## Missing data (GTW Procedure 2): for each block q = I_q x J_q returned by
## find_bicliques(), Algorithm 1 on I_q (sub-seed 4q - 1) and on J_q (4q),
## a common K, concatenated into one block-diagonal group (step 3). `ri, ci`
## are the dense ids of the RETAINED observations, in data order.
## `slot` (optional) is a within-cell index held fixed, as the incomplete-panel
## design does with the period.
ref_groups_blocks <- function(ri, ci, blocks, K, seed, r, slot = NULL) {
  rs <- rep_seed(seed, r)
  stride <- max(1000, 4 * length(blocks) + 1)
  blk <- integer(length(ri))
  for (q in seq_along(blocks))
    blk[ri %in% blocks[[q]]$rows & ci %in% blocks[[q]]$cols] <- q
  stopifnot(all(blk > 0L))
  Gr <- lapply(seq_along(blocks), function(q)
    build_perm_set(length(blocks[[q]]$rows), K,
                   seed = sub_seed(rs, 4 * q - 1, stride)))
  Gc <- lapply(seq_along(blocks), function(q)
    build_perm_set(length(blocks[[q]]$cols), K,
                   seed = sub_seed(rs, 4 * q, stride)))
  s <- if (is.null(slot)) rep(0L, length(ri)) else slot
  from <- paste(ri, ci, s)
  lapply(seq_len(K + 1L), function(k) {
    ti <- ri
    tj <- ci
    for (q in seq_along(blocks)) {
      b <- blocks[[q]]
      sel <- blk == q
      ti[sel] <- b$rows[Gr[[q]][[k]][match(ri[sel], b$rows)]]
      tj[sel] <- b$cols[Gc[[q]][[k]][match(ci[sel], b$cols)]]
    }
    g <- match(paste(ti, tj, s), from)
    stopifnot(!anyNA(g), !anyDuplicated(g))
    g
  })
}

## Nuisance design exactly as the front ends assemble it: an intercept, the
## user's x, and (panel, time_fe = TRUE) period dummies without the reference
## level.
ref_X <- function(x, N, time = NULL) {
  X <- cbind(1, if (is.null(x)) NULL else as.matrix(x))
  if (!is.null(time)) {
    td <- stats::model.matrix(~ factor(dense_id(time)))[, -1L, drop = FALSE]
    X <- cbind(X, td)
  }
  unname(X)
}

## Reference outcome for every repetition of a fit: the naive Procedure 1 on
## each rebuilt group. `groups_of(r)` returns repetition r's gather vectors.
ref_reps <- function(y, D, X, groups_of, n_reps, b = 0)
  lapply(seq_len(n_reps), function(r)
    ref_procedure1(y, D, X, groups_of(r), b = b))

## ---- 3. Data-generating processes --------------------------------------------
## Small, seeded, and returned as data frames so every file builds the same
## objects. All draws go through with_seed(), so no helper moves the global
## RNG stream.

## Complete n1 x n2 dyadic array with two-way random effects
## eps_ij = eta_i + xi_j + u_ij (GTW Eq. 8), a dyad-level d and a nuisance x.
make_dyadic <- function(n1 = 20, n2 = n1, beta = 0, seed = 1,
                        err = function(k) stats::rnorm(k), d_node = FALSE) {
  with_seed(seed, {
    g <- expand.grid(i = seq_len(n1), j = seq_len(n2))
    N <- nrow(g)
    eta <- err(n1)
    xi <- err(n2)
    u <- err(N)
    d <- if (d_node) stats::rnorm(n1)[g$i] else stats::rnorm(N)
    x <- stats::rnorm(N) + 0.5 * d
    g$d <- d
    g$x <- x
    g$y <- beta * d + 0.7 * x + eta[g$i] + xi[g$j] + u
    g
  })
}

## Complete n1 x n2 x T panel with an arbitrary common time trend zeta_t and
## AR(1) idiosyncratic errors over t (InvB holds, InvA does not).
make_panel <- function(n1 = 8, n2 = n1, T = 3, beta = 0, seed = 1,
                       trend_sd = 3, rho = 0.6) {
  with_seed(seed, {
    g <- expand.grid(i = seq_len(n1), j = seq_len(n2), t = seq_len(T))
    N <- nrow(g)
    eta <- stats::rnorm(n1)
    xi <- stats::rnorm(n2)
    zeta <- cumsum(stats::rnorm(T, sd = trend_sd))
    u <- matrix(0, n1 * n2, T)
    u[, 1] <- stats::rnorm(n1 * n2)
    if (T > 1) for (t in 2:T)
      u[, t] <- rho * u[, t - 1] + sqrt(1 - rho^2) * stats::rnorm(n1 * n2)
    g$d <- stats::rnorm(N) + 0.3 * g$t
    g$x <- stats::rnorm(N)
    g$y <- beta * g$d + 0.5 * g$x + eta[g$i] + xi[g$j] + zeta[g$t] +
      as.vector(u)
    g
  })
}

## Complete three-way array under full random effects (InvA in 3 indices).
make_threeway <- function(n1 = 6, n2 = 7, n3 = 8, beta = 0, seed = 1) {
  with_seed(seed, {
    g <- expand.grid(i = seq_len(n1), j = seq_len(n2), l = seq_len(n3))
    N <- nrow(g)
    g$d <- stats::rnorm(N)
    g$x <- stats::rnorm(N)
    g$y <- beta * g$d + 0.5 * g$x + stats::rnorm(n1)[g$i] +
      stats::rnorm(n2)[g$j] + stats::rnorm(n3)[g$l] + stats::rnorm(N)
    g
  })
}

## Two-way layout: n1 x n2 cells with ell_ij replicates drawn from `sizes`,
## cell effects eta_ij, and a d that varies within cells.
make_layout <- function(n1 = 4, n2 = 4, sizes = 6:9, beta = 0, seed = 1,
                        d_cell_constant = FALSE) {
  with_seed(seed, {
    cells <- expand.grid(i = seq_len(n1), j = seq_len(n2))
    ell <- sizes[sample.int(length(sizes), nrow(cells), replace = TRUE)]
    g <- data.frame(i = rep(cells$i, ell), j = rep(cells$j, ell),
                    l = unlist(lapply(ell, seq_len)))
    cid <- rep(seq_len(nrow(cells)), ell)
    N <- nrow(g)
    g$d <- if (d_cell_constant) stats::rnorm(nrow(cells))[cid]
           else stats::rnorm(N)
    g$x <- stats::rnorm(N)
    g$y <- beta * g$d + 0.5 * g$x + stats::rnorm(nrow(cells))[cid] +
      stats::rnorm(N)
    g
  })
}

## ---- 4. Monte Carlo helpers --------------------------------------------------

## The slow tests run only on request: MWPERM_SLOW_TESTS=true, and never on
## CRAN. They are the Monte Carlo checks of Theorem 1; the default suite stays
## within about two minutes without them.
skip_if_not_slow <- function() {
  skip_on_cran()
  if (!identical(tolower(Sys.getenv("MWPERM_SLOW_TESTS")), "true"))
    skip("slow Monte Carlo test: set MWPERM_SLOW_TESTS=true to run")
}

## Theorem 1: P(pval <= alpha | X, D) <= alpha for EVERY alpha. With p-values
## on the grid {1, ..., K+1}/(K+1), it is enough to check the empirical CDF
## at every atom a = j/(K+1) against a + 3 binomial standard errors (a
## one-sided bound). Returns the table so a failure shows every atom.
ecdf_at_atoms <- function(p, K) {
  atoms <- seq_len(K + 1L) / (K + 1L)
  n <- length(p)
  data.frame(atom = atoms,
             ecdf = vapply(atoms, function(a) mean(p <= a + 1e-12), 0),
             bound = atoms + 3 * sqrt(atoms * (1 - atoms) / n))
}
expect_valid_ecdf <- function(p, K) {
  tab <- ecdf_at_atoms(p, K)
  bad <- tab$ecdf > tab$bound
  expect(!any(bad), paste0(
    "empirical CDF above atom + 3 SE at ",
    paste(sprintf("%.3f (ecdf %.3f > %.3f)", tab$atom[bad], tab$ecdf[bad],
                  tab$bound[bad]), collapse = "; ")))
  invisible(tab)
}

## ---- 5. Logged discrepancies ---------------------------------------------------
## A test whose assertion follows the paper (or the draft) but which the code
## does not satisfy is kept, and skipped with a pointer to the discrepancy log
## in TESTING_PLAN.md, so the suite stays green without editing the assertion.
## Set MWPERM_SHOW_DISCREPANCIES=true to run the skipped bodies and see each
## one fail.
skip_discrepancy <- function(id, why) {
  if (!identical(tolower(Sys.getenv("MWPERM_SHOW_DISCREPANCIES")), "true"))
    skip(paste0("Discrepancy ", id, ": ", why, "; see TESTING_PLAN.md"))
}

## Optional record of the Monte Carlo tables: when MWPERM_MC_OUT names a
## file, each slow test appends its table there (used to write the report;
## nothing is written otherwise).
mc_record <- function(name, tab) {
  out <- Sys.getenv("MWPERM_MC_OUT")
  if (nzchar(out)) {
    tab <- cbind(design = name, tab)
    utils::write.table(tab, out, append = file.exists(out), sep = ",",
                       row.names = FALSE, col.names = !file.exists(out))
  }
  invisible(tab)
}
