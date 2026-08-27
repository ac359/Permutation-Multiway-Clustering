## Every permutation the engine sees must be a bijection of the observations.
##
## The permutation builders key an observation by its cluster coordinates and
## translate a permuted cell code back to a row index with match() (or a
## position table). Both resolve a repeated cell to its FIRST observation, so a
## duplicated cell turns the gather vector into a many-to-one map: the statistic
## would then be computed on duplicated rows, with no NA and no warning
## anywhere. These assertions convert that silent wrong answer into an error.
##
## Reaches internals via ::: -- run against a FRESHLY INSTALLED package.
library(mwperm)

ok <- function(...) stopifnot(...)
err <- function(expr, pat) {
  m <- tryCatch({ expr; NA_character_ },
                error = function(e) conditionMessage(e))
  ok(!is.na(m), grepl(pat, m, fixed = TRUE))
  invisible(m)
}

## ---- 1. duplicated cells rejected at the builder entry ---------------------
G <- build_perm_set(3L, 2L, seed = 1)
dup_coords <- cbind(c(1L, 1L, 2L, 2L, 3L, 3L), c(1L, 1L, 2L, 2L, 3L, 3L))
err(mwperm:::.build_obs_perms(dup_coords, list(G, G), design = "dyadic"),
    "requires exactly one observation per cell")
## the message names the design and points somewhere useful
m <- tryCatch(mwperm:::.build_obs_perms(dup_coords, list(G, G),
                                        design = "dyadic",
                                        front_end = "mwperm_layout()"),
              error = function(e) conditionMessage(e))
ok(grepl("dyadic", m, fixed = TRUE), grepl("mwperm_layout()", m, fixed = TRUE))

## the block builder guards the same invariant
err(mwperm:::.build_obs_perms_blocks(
      1L, 2L, list(list(rows = 1:3, cols = 1:3)),
      ri = c(1L, 1L, 2L, 3L, 2L, 3L), ci = c(1L, 1L, 2L, 2L, 3L, 3L),
      blk = rep(1L, 6L), lrow = c(1L, 1L, 2L, 3L, 2L, 3L),
      lcol = c(1L, 1L, 2L, 2L, 3L, 3L)),
    "exactly one observation per (row, col) cell")

## ---- 2. the bijection assertion itself --------------------------------------
err(mwperm:::.assert_bijection(list(c(1L, 2L, 2L)), 3L, "dyadic"),
    "is not a bijection")
err(mwperm:::.assert_bijection(list(1:2), 3L, "panel"), "is not a bijection")
err(mwperm:::.assert_bijection(list(c(1L, NA_integer_, 3L)), 3L, "threeway"),
    "is not a bijection")
ok(is.null(mwperm:::.assert_bijection(list(c(3L, 1L, 2L), 1:3), 3L, "ok")))

## ---- 3. every builder in the package returns bijections --------------------
## dyadic / panel / threeway, via .build_obs_perms
gg <- expand.grid(i = 1:5, j = 1:5)
Gd <- build_perm_set(5L, 4L, seed = 2)
for (grp in list(list(Gd, Gd), list(Gd, NULL), list(NULL, Gd))) {
  ops <- mwperm:::.build_obs_perms(cbind(gg$i, gg$j), grp)
  for (o in ops) ok(length(o) == nrow(gg), !anyDuplicated(o),
                    identical(sort(o), seq_len(nrow(gg))))
}
## panel: three coordinates, time held fixed
tp <- expand.grid(i = 1:4, j = 1:4, t = 1:3)
G4 <- build_perm_set(4L, 3L, seed = 3)
ops <- mwperm:::.build_obs_perms(cbind(tp$i, tp$j, tp$t), list(G4, G4, NULL))
for (o in ops) {
  ok(identical(sort(o), seq_len(nrow(tp))))
  ok(identical(tp$t[o], tp$t))          # time really is held fixed
}

## layout: within-cell permutations
cellv <- rep(1:6, each = 4L)
widx <- rep(1:4, times = 6L)
cg <- lapply(1:6, function(k) build_perm_set(4L, 3L, seed = 100 + k))
ops <- mwperm:::.build_obs_perms_layout(cellv, widx, cg)
for (o in ops) {
  ok(identical(sort(o), seq_along(cellv)))
  ok(identical(cellv[o], cellv))        # permutation stays inside its cell
}

## missing / irregular: block-diagonal, with and without slots
blocks <- list(list(rows = 1:3, cols = 1:3), list(rows = 4:6, cols = 4:6))
ri <- c(rep(1:3, each = 3L), rep(4:6, each = 3L))
ci <- c(rep(1:3, 3L), rep(4:6, 3L))
blk <- rep(1:2, each = 9L)
lrow <- c(rep(1:3, each = 3L), rep(1:3, each = 3L))
lcol <- c(rep(1:3, 3L), rep(1:3, 3L))
ops <- mwperm:::.build_obs_perms_blocks(5L, 2L, blocks, ri, ci, blk, lrow, lcol)
for (o in ops) ok(identical(sort(o), seq_along(ri)))

## with a within-cell slot: the slot must be held fixed and the SAME (pi, sigma)
## must act in every slot
ri2 <- rep(ri, 2L); ci2 <- rep(ci, 2L); blk2 <- rep(blk, 2L)
lrow2 <- rep(lrow, 2L); lcol2 <- rep(lcol, 2L)
slot <- rep(1:2, each = length(ri))
ops <- mwperm:::.build_obs_perms_blocks(5L, 2L, blocks, ri2, ci2, blk2,
                                        lrow2, lcol2, slot = slot)
for (o in ops) {
  ok(identical(sort(o), seq_along(ri2)))
  ok(identical(slot[o], slot))
  ok(identical(cbind(ri2[o], ci2[o])[slot == 1L, ],
               cbind(ri2[o], ci2[o])[slot == 2L, ]))
}
## a repeated (row, col, slot) key is an internal error, not a silent answer
err(mwperm:::.build_obs_perms_blocks(5L, 2L, blocks, ri2, ci2, blk2,
                                     lrow2, lcol2, slot = rep(1L, length(ri2))),
    "appears more than once")

## ---- 4. front ends still reject repeated cells with a usable message -------
set.seed(6)
dd <- expand.grid(i = 1:4, j = 1:4)
dd <- rbind(dd, dd[1L, ])                      # one duplicated cell
err(mwperm_dyadic(rnorm(nrow(dd)), rnorm(nrow(dd)), row = dd$i, col = dd$j,
                  conf_int = FALSE, seed = 1),
    "one observation per (row, col) cell")

cat("test-bijection.R: all assertions passed\n")
