## Phase 2.2 — the OBSERVATION-level permutation set must itself be a cyclic group
## of order K+1 for every design. Closure of this joint group is what Theorem 1
## rests on (GTW 2026, Proposition 2 lifted to observations); if it fails for a
## design the finite-sample guarantee is void there. Also: a permutation to an
## unobserved cell must ERROR, never silently NA.
library(mwperm)

compose <- function(a, b) a[b]                 # gather-vector composition
## O is closed as a cyclic group iff O[[r]] o O[[s]] == O[[(r+s) mod B]] for all r,s
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

## ---- 1. dyadic: paired row/col groups compose cyclically ---------------------
g <- expand.grid(i = 1:6, j = 1:6); co <- cbind(g$i, g$j)
Od <- bop(co, list(build_perm_set(6, K, seed = 1), build_perm_set(6, K, seed = 2)))
stopifnot(is_closed(Od))
stopifnot(identical(Od[[1]], seq_len(nrow(co))))          # element 0 = identity

## ---- 2. three-way (InvA): all three dimensions permuted jointly --------------
g3 <- expand.grid(i = 1:5, j = 1:5, l = 1:5); co3 <- cbind(g3$i, g3$j, g3$l)
O3 <- bop(co3, list(build_perm_set(5, K, seed = 1), build_perm_set(5, K, seed = 2),
                    build_perm_set(5, K, seed = 3)))
stopifnot(is_closed(O3))

## ---- 3. panel (InvB): SAME (pi,sigma) every period, time held fixed ----------
gp <- expand.grid(i = 1:5, j = 1:5, t = 1:4); cop <- cbind(gp$i, gp$j, gp$t)
Op <- bop(cop, list(build_perm_set(5, K, seed = 1), build_perm_set(5, K, seed = 2), NULL))
stopifnot(is_closed(Op))
## the third (time) coordinate is untouched by every element -> any time trend
## is preserved and partialled out, exactly as InvB requires
stopifnot(all(vapply(Op, function(o) identical(cop[o, 3], cop[, 3]), logical(1))))

## ---- 4. two-way layout: within-cell permutations only ------------------------
gl <- expand.grid(rep = 1:4, i = 1:3, j = 1:3)
cell <- as.integer(interaction(gl$i, gl$j, drop = TRUE)); ncell <- max(cell)
widx <- integer(length(cell))
for (c in seq_len(ncell)) {
  idx <- which(cell == c); widx[idx] <- rank(gl$rep[idx], ties.method = "first")
}
cg <- lapply(seq_len(ncell), function(c) build_perm_set(4L, K, seed = 100 + c))
Ol <- bopl(cell, widx, cg)
stopifnot(is_closed(Ol))
## no observation is ever moved out of its (i,j) cell
stopifnot(all(vapply(Ol, function(o) identical(cell[o], cell), logical(1))))

## ---- 5. missing block-diagonal: two disjoint blocks --------------------------
## The block-diagonal builder is defined inline in mwperm_missing() and is not
## independently callable (see FINDING F2.6 in audit/02_correctness.md). We
## reconstruct it verbatim to pin its group closure; the reconstruction mirrors
## R/missing.R:131-164 and shares that logic (this is a closure check, not the
## independent oracle of 2.3).
blocks <- list(list(rows = 1:4, cols = 1:4), list(rows = 5:8, cols = 5:8))
gg <- rbind(expand.grid(r = 1:4, c = 1:4), expand.grid(r = 5:8, c = 5:8))
ri_k <- gg$r; ci_k <- gg$c; blk_k <- ifelse(ri_k <= 4, 1L, 2L)
lrow_k <- integer(nrow(gg)); lcol_k <- integer(nrow(gg))
for (q in 1:2) {
  sel <- blk_k == q
  lrow_k[sel] <- match(ri_k[sel], blocks[[q]]$rows)
  lcol_k[sel] <- match(ci_k[sel], blocks[[q]]$cols)
}
code_k <- cc(cbind(ri_k, ci_k))
build_missing_ops <- function(seed) {
  rowG <- colG <- vector("list", 2)
  for (q in 1:2) {
    rowG[[q]] <- build_perm_set(4L, K, seed = seed * 1000L + 4L * q - 1L)
    colG[[q]] <- build_perm_set(4L, K, seed = seed * 1000L + 4L * q)
  }
  ops <- vector("list", K + 1L)
  for (k in seq_len(K + 1L)) {
    tg_row <- ri_k; tg_col <- ci_k
    for (q in 1:2) {
      sel <- blk_k == q
      tg_row[sel] <- blocks[[q]]$rows[rowG[[q]][[k]][lrow_k[sel]]]
      tg_col[sel] <- blocks[[q]]$cols[colG[[q]][[k]][lcol_k[sel]]]
    }
    ops[[k]] <- match(cc(cbind(tg_row, tg_col)), code_k)
  }
  ops
}
Om <- build_missing_ops(1L)
stopifnot(is_closed(Om))
stopifnot(!anyNA(unlist(Om)))                              # never leaves the observed set
stopifnot(all(vapply(Om, function(o) identical(blk_k[o], blk_k), logical(1))))  # stays in-block

## ---- 6. a permutation to an unobserved cell ERRORS (not NA) ------------------
gi <- expand.grid(i = 1:5, j = 1:5); coi <- cbind(gi$i, gi$j)[-1, ]   # drop cell (1,1)
msg <- tryCatch({
  bop(coi, list(build_perm_set(5, K, seed = 1), build_perm_set(5, K, seed = 2)))
  NA_character_
}, error = function(e) conditionMessage(e))
stopifnot(!is.na(msg), grepl("unobserved", msg))

cat("test-obsperms.R: all assertions passed\n")
