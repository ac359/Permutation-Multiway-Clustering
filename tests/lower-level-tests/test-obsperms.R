## The OBSERVATION-level permutation set must itself be a cyclic group
## of order K+1 for every design. Closure of this joint group is what Theorem 1
## rests on (GTW 2026, Proposition 2 lifted to observations); if it fails for a
## design the finite-sample guarantee is void there. Also: a permutation to an
## unobserved cell must ERROR, never silently NA.
##
## Lower-level test (tests/lower-level-tests/): it exercises the machinery
## beneath the front ends, mostly through package internals, so it must run
## against a FRESHLY INSTALLED copy -- mwperm::: resolves against the
## installed package, never against a source()d working tree.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

compose <- function(a, b) a[b]                 # gather-vector composition
## O is closed as a cyclic group iff O[[r]] o O[[s]] == O[[(r+s) mod B]] for all
## r,s
is_closed <- function(O) {
  B <- length(O)
  for (r in 0:(B - 1L)) for (s in 0:(B - 1L))
    if (!identical(compose(O[[r + 1L]], O[[s + 1L]]), O[[((r + s) %% B) + 1L]]))
      return(FALSE)
  TRUE
}
K <- 3L
bop  <- mwperm:::.build_obs_perms
bopl <- mwperm:::.build_obs_perms_layout
cc   <- mwperm:::.cell_code

## ---- 1. dyadic: paired row/col groups compose cyclically ------------------
g <- expand.grid(i = 1:6, j = 1:6)
co <- cbind(g$i, g$j)
Od <- bop(co, list(build_perm_set(6, K, seed = 1), build_perm_set(6, K,
                                                                  seed = 2)))
stopifnot(is_closed(Od))
stopifnot(identical(Od[[1]], seq_len(nrow(co))))          # element 0 = identity

## ---- 2. three-way (InvA): all three dimensions permuted jointly -----------
g3 <- expand.grid(i = 1:5, j = 1:5, l = 1:5)
co3 <- cbind(g3$i, g3$j, g3$l)
O3 <- bop(co3, list(build_perm_set(5, K, seed = 1), build_perm_set(5, K,
                                                                   seed = 2),
                    build_perm_set(5, K, seed = 3)))
stopifnot(is_closed(O3))

## ---- 3. panel (InvB): SAME (pi,sigma) every period, time held fixed -------
gp <- expand.grid(i = 1:5, j = 1:5, t = 1:4)
cop <- cbind(gp$i, gp$j, gp$t)
Op <- bop(cop, list(build_perm_set(5, K, seed = 1),
                    build_perm_set(5, K, seed = 2), NULL))
stopifnot(is_closed(Op))
## the third (time) coordinate is untouched by every element -> any time trend
## is preserved and partialled out, exactly as InvB requires
stopifnot(all(vapply(Op, function(o) identical(cop[o, 3], cop[, 3]),
                     logical(1))))

## ---- 4. two-way layout: within-cell permutations only ---------------------
gl <- expand.grid(rep = 1:4, i = 1:3, j = 1:3)
cell <- as.integer(interaction(gl$i, gl$j, drop = TRUE))
ncell <- max(cell)
widx <- integer(length(cell))
for (c in seq_len(ncell)) {
  idx <- which(cell == c)
  widx[idx] <- rank(gl$rep[idx], ties.method = "first")
}
cg <- lapply(seq_len(ncell), function(c) build_perm_set(4L, K, seed = 100 + c))
Ol <- bopl(cell, widx, cg)
stopifnot(is_closed(Ol))
## no observation is ever moved out of its (i,j) cell
stopifnot(all(vapply(Ol, function(o) identical(cell[o], cell), logical(1))))

## ---- 5. missing block-diagonal: two disjoint blocks -----------------------
## Closure of the SAME code the fit runs: the builder was extracted from
## mwperm_missing()'s inline closure to the named internal
## .build_obs_perms_blocks(), so this test can no longer
## drift from the production code path.
blocks <- list(list(rows = 1:4, cols = 1:4), list(rows = 5:8, cols = 5:8))
gg <- rbind(expand.grid(r = 1:4, c = 1:4), expand.grid(r = 5:8, c = 5:8))
ri_k <- gg$r
ci_k <- gg$c
blk_k <- ifelse(ri_k <= 4, 1L, 2L)
lrow_k <- integer(nrow(gg))
lcol_k <- integer(nrow(gg))
for (q in 1:2) {
  sel <- blk_k == q
  lrow_k[sel] <- match(ri_k[sel], blocks[[q]]$rows)
  lcol_k[sel] <- match(ci_k[sel], blocks[[q]]$cols)
}
Om <- mwperm:::.build_obs_perms_blocks(1L, K, blocks, ri = ri_k, ci = ci_k,
                                       blk = blk_k, lrow = lrow_k,
                                       lcol = lcol_k)
stopifnot(is_closed(Om))
stopifnot(!anyNA(unlist(Om)))                # never leaves the observed set
stopifnot(all(vapply(Om, function(o) identical(blk_k[o], blk_k),
                     logical(1))))  # stays in-block

## ---- 6. a permutation to an unobserved cell ERRORS (not NA) ---------------
gi <- expand.grid(i = 1:5, j = 1:5)
coi <- cbind(gi$i, gi$j)[-1, ]   # drop cell (1,1)
msg <- tryCatch({
  bop(coi, list(build_perm_set(5, K, seed = 1), build_perm_set(5, K, seed = 2)))
  NA_character_
}, error = function(e) conditionMessage(e))
stopifnot(!is.na(msg), grepl("unobserved", msg))

## ---- 7. the position-table and match() branches are the same map ---------
## Both builders translate permuted cell codes back to observation indices
## either through a position table over the whole mixed-radix index space
## (fast, but that space can be enormous for sparse ids) or through match().
## The branch is chosen by size alone and must never change a gather vector;
## `pos_table` forces each branch so the equivalence is asserted, not assumed.
## Dyadic, dense ids: both branches, and the automatic choice, agree.
grp <- list(build_perm_set(6, K, seed = 1), build_perm_set(6, K, seed = 2))
stopifnot(identical(bop(co, grp, pos_table = TRUE),
                    bop(co, grp, pos_table = FALSE)),
          identical(bop(co, grp), Od))
## Panel keying (a held-fixed third digit) through both branches too.
gp <- expand.grid(i = 1:5, j = 1:5, t = 1:3)
cop <- cbind(gp$i, gp$j, gp$t)
gsp <- list(build_perm_set(5, K, seed = 3), build_perm_set(5, K, seed = 4),
            NULL)
stopifnot(identical(bop(cop, gsp, pos_table = TRUE),
                    bop(cop, gsp, pos_table = FALSE)),
          identical(bop(cop, gsp), bop(cop, gsp, pos_table = TRUE)))
## The gate itself: a complete array always has n_cells == N, so the table is
## taken unless the array has more than 2^24 cells; a block design over sparse
## global ids can have an index space far larger than the retained data, and
## then match() is used.
upt <- mwperm:::.use_pos_table
stopifnot(upt(1600, 1600), upt(2^24, 2^24), !upt(2^24 + 1, 2^24 + 1),
          upt(64 * 100, 100), !upt(64 * 100 + 1, 100))
## The block builder: same contract, both branches forced.
stopifnot(identical(mwperm:::.build_obs_perms_blocks(
                      1L, K, blocks, ri = ri_k, ci = ci_k, blk = blk_k,
                      lrow = lrow_k, lcol = lcol_k, pos_table = TRUE),
                    mwperm:::.build_obs_perms_blocks(
                      1L, K, blocks, ri = ri_k, ci = ci_k, blk = blk_k,
                      lrow = lrow_k, lcol = lcol_k, pos_table = FALSE)),
          identical(mwperm:::.build_obs_perms_blocks(
                      1L, K, blocks, ri = ri_k, ci = ci_k, blk = blk_k,
                      lrow = lrow_k, lcol = lcol_k), Om))
## and with a slot (the irregular / incomplete-panel keying)
sl <- rep(1:2, each = nrow(gg))
stopifnot(identical(mwperm:::.build_obs_perms_blocks(
                      1L, K, blocks, ri = rep(ri_k, 2), ci = rep(ci_k, 2),
                      blk = rep(blk_k, 2), lrow = rep(lrow_k, 2),
                      lcol = rep(lcol_k, 2), slot = sl, pos_table = TRUE),
                    mwperm:::.build_obs_perms_blocks(
                      1L, K, blocks, ri = rep(ri_k, 2), ci = rep(ci_k, 2),
                      blk = rep(blk_k, 2), lrow = rep(lrow_k, 2),
                      lcol = rep(lcol_k, 2), slot = sl, pos_table = FALSE)))
## Sparse global ids: the same two 4 x 4 blocks placed at ids near 1 and near
## 5000, so the index space is 25 million cells for 32 observations. The
## automatic choice must be match(), and its output must equal the forced
## table branch on the DENSE relabelling of the same design (the block-local
## structure, and hence the gather vectors, are identical by construction).
big <- c(1:4, 4997:5000)
blocks_sp <- list(list(rows = big[1:4], cols = big[1:4]),
                  list(rows = big[5:8], cols = big[5:8]))
O_sp <- mwperm:::.build_obs_perms_blocks(1L, K, blocks_sp, ri = big[ri_k],
                                         ci = big[ci_k], blk = blk_k,
                                         lrow = lrow_k, lcol = lcol_k)
stopifnot(!upt(5000 * 5000, nrow(gg)),
          identical(O_sp, mwperm:::.build_obs_perms_blocks(
                      1L, K, blocks_sp, ri = big[ri_k], ci = big[ci_k],
                      blk = blk_k, lrow = lrow_k, lcol = lcol_k,
                      pos_table = FALSE)),
          identical(O_sp, Om))

passed("test-obsperms.R")
