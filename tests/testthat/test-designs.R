## The designs of Guo, Toulis & Wang (2026), Sections 5-6 (draft Section
## 2.7): the same engine run with different groups. For each design this file
## checks (a) the STRUCTURE of the group the package builds -- which
## coordinates move, and how -- on the package's own gather vectors, and
## (b) that every per-repetition p-value equals the naive Procedure 1 of
## helper-reference.R run with the group the paper prescribes, rebuilt from
## the seed. (b) is what would catch a structurally wrong group that still
## happened to be a group.

expect_reps_match <- function(fit, y, D, X, groups_of) {
  ref <- ref_reps(y, D, X, groups_of, length(fit$pvalues_rep))
  for (r in seq_along(ref)) expect_gt(ref[[r]]$gap, 1e-6)
  expect_identical(fit$pvalues_rep, vapply(ref, `[[`, 0, "pvalue"))
}

## ---- three-way (Section 6.1, InvA) -------------------------------------------

test_that("three-way: independent groups per dimension, applied jointly", {
  tw <- make_threeway(6, 7, 8, beta = 0.15, seed = 81)
  ids <- list(tw$i, tw$j, tw$l)
  fit <- mwperm_threeway(tw$y, tw$d, x = tw$x, id1 = tw$i, id2 = tw$j,
                         id3 = tw$l, n_reps = 3, seed = 9, conf_int = FALSE)
  expect_identical(fit$K, 5L)                       # min(6, 7, 8) - 1
  expect_reps_match(fit, tw$y, tw$d, ref_X(tw$x, nrow(tw)),
                    function(r) ref_groups_crossed(ids, 5, 9, r))
  ## structure: every non-identity element moves every coordinate, by three
  ## different relabellings (sub-seeds 1, 2, 3 of the repetition seed)
  rs <- rep_seed(9, 1)
  G <- lapply(1:3, function(j) build_perm_set(c(6, 7, 8)[j], 5,
                                              seed = sub_seed(rs, j)))
  pkg <- mwperm:::.build_obs_perms(cbind(tw$i, tw$j, tw$l), G)
  for (k in 2:6) {
    g <- pkg[[k]]
    expect_identical(tw$i[g], G[[1]][[k]][tw$i])
    expect_identical(tw$j[g], G[[2]][[k]][tw$j])
    expect_identical(tw$l[g], G[[3]][[k]][tw$l])
  }
  expect_false(identical(G[[1]][[2]][1:6], G[[2]][[2]][1:6]))
})

## ---- panel (Section 6.2, InvB) -----------------------------------------------

test_that("panel: ONE (pi, sigma) shared by every period, time held fixed", {
  pn <- make_panel(9, 8, T = 4, beta = 0.1, seed = 82)
  rs <- rep_seed(3, 1)
  Gr <- build_perm_set(9, 7, seed = sub_seed(rs, 1))
  Gc <- build_perm_set(8, 7, seed = sub_seed(rs, 2))
  pkg <- mwperm:::.build_obs_perms(cbind(pn$i, pn$j, pn$t),
                                   list(Gr, Gc, NULL))
  for (k in 2:8) {
    g <- pkg[[k]]
    expect_identical(pn$t[g], pn$t)                         # period fixed
    expect_identical(pn$i[g], Gr[[k]][pn$i])                # same pi ...
    expect_identical(pn$j[g], Gc[[k]][pn$j])                # ... and sigma
  }
})

test_that("panel: p-values match Procedure 1 with period dummies in X", {
  ## time_fe = TRUE adds the period dummies to X; they are invariant under
  ## the within-period permutation, so [X | X_k] is rank deficient -- the
  ## reference handles it by numerical rank.
  pn <- make_panel(9, 8, T = 4, beta = 0.1, seed = 82)
  ids <- list(pn$i, pn$j, pn$t)
  groups <- function(r) ref_groups_crossed(ids, 7, 3, r, permuted = 1:2)
  X <- ref_X(pn$x, nrow(pn), time = pn$t)
  M <- cbind(X, X[groups(1)[[2]], ])
  expect_lt(qr(M)$rank, ncol(M) - 1L)     # intercept AND dummies duplicated
  fit <- mwperm_panel(pn$y, pn$d, x = pn$x, row = pn$i, col = pn$j,
                      time = pn$t, n_reps = 3, seed = 3, conf_int = FALSE)
  expect_reps_match(fit, pn$y, pn$d, X, groups)
  ## and without the dummies
  fit0 <- mwperm_panel(pn$y, pn$d, x = pn$x, row = pn$i, col = pn$j,
                       time = pn$t, n_reps = 3, seed = 3, conf_int = FALSE,
                       time_fe = FALSE)
  expect_reps_match(fit0, pn$y, pn$d, ref_X(pn$x, nrow(pn)), groups)
})

## ---- two-way layout (Section 6.3) --------------------------------------------

test_that("layout: permutations act within cells, over the replicate index", {
  ly <- make_layout(4, 4, sizes = 6:9, beta = 0.3, seed = 83)
  fit <- mwperm_layout(ly$y, ly$d, x = ly$x, row = ly$i, col = ly$j,
                       rep = ly$l, n_reps = 3, seed = 4, conf_int = FALSE)
  expect_identical(fit$K, 5L)                        # smallest cell - 1
  groups <- function(r) ref_groups_layout(ly$i, ly$j, ly$l, 5, 4, r)
  expect_reps_match(fit, ly$y, ly$d, ref_X(ly$x, nrow(ly)), groups)
  ## structure, on the package's own gather vectors: nothing leaves its cell,
  ## and inside cell c the map is that cell's Algorithm 1 group
  cell <- dense_id(interaction(dense_id(ly$i), dense_id(ly$j), drop = TRUE))
  slot <- ref_slot(cell, ly$l)
  rs <- rep_seed(4, 1)
  ell <- tabulate(cell)
  G <- lapply(seq_along(ell), function(c)
    build_perm_set(ell[c], 5, seed = sub_seed(rs, c)))
  pkg <- mwperm:::.build_obs_perms_layout(cell, slot, G)
  for (k in 2:6) {
    g <- pkg[[k]]
    expect_identical(cell[g], cell)
    expect_identical(slot[g], vapply(seq_along(cell), function(o)
      G[[cell[o]]][[k]][slot[o]], 1L))
  }
})

test_that("layout L0: keeps cells with >= L0 replicates, each cut to L0", {
  ly <- make_layout(5, 5, sizes = c(3, 5, 8, 11), beta = 0.3, seed = 84)
  cell <- paste(ly$i, ly$j)
  ell <- table(cell)
  L0 <- 5
  fit <- mwperm_layout(ly$y, ly$d, x = ly$x, row = ly$i, col = ly$j,
                       rep = ly$l, L0 = L0, n_reps = 2, seed = 6,
                       conf_int = FALSE)
  expect_identical(fit$n_obs, as.integer(L0 * sum(ell >= L0)))
  expect_identical(fit$K, as.integer(L0 - 1))
  expect_identical(unname(fit$n_clusters["ncell"]), as.integer(sum(ell >= L0)))
  expect_match(fit$note[1], sprintf("kept %d of %d cells", sum(ell >= L0),
                                    length(ell)), fixed = TRUE)
  ## the retained rows: the draw .downsample_to_L0() makes from `seed`
  keep <- mwperm:::.downsample_to_L0(dense_id(interaction(
    dense_id(ly$i), dense_id(ly$j), drop = TRUE)),
    as.vector(tabulate(dense_id(interaction(dense_id(ly$i), dense_id(ly$j),
                                            drop = TRUE)))), L0, seed = 6)
  kept <- ly[keep, ]
  expect_true(all(table(paste(kept$i, kept$j)) == L0))
  groups <- function(r) ref_groups_layout(kept$i, kept$j, kept$l, L0 - 1, 6, r)
  expect_reps_match(fit, kept$y, kept$d, ref_X(kept$x, nrow(kept)), groups)
})

test_that("layout under a shared replicate effect zeta_l (GTW Sec. 6.3)", {
  ## GTW Section 6.3 (and the draft's Section 2.7) state the layout test is
  ## valid for eps_ijl = eta_ij + zeta_l + u_ijl with zeta_l i.i.d. and SHARED
  ## by all cells. For the group to leave that law invariant, element k must
  ## apply the same permutation of l in every cell (otherwise the cross-cell
  ## alignment of zeta changes). Algorithm 1 run independently per cell, as
  ## step (i) prescribes and the package does, does not do that.
  skip_discrepancy("D7", paste("GTW Sec. 6.3 claims validity under a shared",
                               "zeta_l, but per-cell Algorithm 1 draws apply",
                               "different permutations of l in different",
                               "cells, so that law is not invariant"))
  ly <- make_layout(3, 3, sizes = 6, seed = 92)       # equal cell sizes
  cell <- dense_id(interaction(dense_id(ly$i), dense_id(ly$j), drop = TRUE))
  rs <- rep_seed(1, 1)
  G <- lapply(1:9, function(c) build_perm_set(6, 5, seed = sub_seed(rs, c)))
  for (k in 2:6) {
    maps <- lapply(G, `[[`, k)
    expect_true(all(vapply(maps, identical, TRUE, maps[[1]])))
  }
})

test_that("layout: a cell-constant regressor draws the no-power warning", {
  ly <- make_layout(4, 4, sizes = 6:8, seed = 85, d_cell_constant = TRUE)
  expect_warning(
    mwperm_layout(ly$y, ly$d, row = ly$i, col = ly$j, rep = ly$l,
                  n_reps = 1, seed = 1, conf_int = FALSE),
    "constant within every cell")
})

## ---- missing data (Section 5, Procedure 2) -----------------------------------

test_that("missing: Procedure 2 on the pooled blocks, K + 1 <= smallest side", {
  dy <- make_dyadic(14, 14, beta = 0.2, seed = 86)
  obs <- with_seed(87, stats::runif(nrow(dy)) < 0.8)       # MCAR mask
  md <- dy[obs, ]
  fit <- mwperm_missing(md$y, md$d, x = md$x, row = md$i, col = md$j,
                        min_block = 3, n_reps = 3, seed = 7, conf_int = FALSE)
  blocks <- find_bicliques(dense_id(md$i), dense_id(md$j), min_block = 3)
  sides <- vapply(blocks, function(b) min(length(b$rows), length(b$cols)), 1L)
  expect_identical(fit$n_blocks, length(blocks))
  expect_identical(fit$K, min(sides) - 1L)
  cells <- sum(vapply(blocks, function(b) length(b$rows) * length(b$cols), 1))
  expect_identical(fit$cells_used, as.integer(cells))
  expect_identical(fit$n_obs, as.integer(cells))
  expect_identical(fit$cells_total, nrow(md))
  expect_identical(nobs(fit), fit$n_obs)
  ## the retained observations, in data order, and the block-diagonal group
  ri <- dense_id(md$i)
  ci <- dense_id(md$j)
  keep <- Reduce(`|`, lapply(blocks, function(b)
    ri %in% b$rows & ci %in% b$cols))
  kd <- md[keep, ]
  groups <- function(r) ref_groups_blocks(ri[keep], ci[keep], blocks,
                                          fit$K, 7, r)
  expect_reps_match(fit, kd$y, kd$d, ref_X(kd$x, nrow(kd)), groups)
  ## a K the smallest block cannot carry is refused
  expect_error(mwperm_missing(md$y, md$d, x = md$x, row = md$i, col = md$j,
                              min_block = 3, K = min(sides), seed = 1),
               "cannot be larger than that")
})

test_that("missing: blocks are permuted separately and never mixed", {
  dy <- make_dyadic(12, 12, seed = 88)
  nd <- subset(dy, i != j)
  ri <- dense_id(nd$i)
  ci <- dense_id(nd$j)
  blocks <- find_bicliques(ri, ci, min_block = 3)
  expect_gte(length(blocks), 2L)
  keep <- Reduce(`|`, lapply(blocks, function(b)
    ri %in% b$rows & ci %in% b$cols))
  blk <- integer(sum(keep))
  for (q in seq_along(blocks))
    blk[ri[keep] %in% blocks[[q]]$rows & ci[keep] %in% blocks[[q]]$cols] <- q
  lrow <- lcol <- integer(sum(keep))
  for (q in seq_along(blocks)) {
    s <- blk == q
    lrow[s] <- match(ri[keep][s], blocks[[q]]$rows)
    lcol[s] <- match(ci[keep][s], blocks[[q]]$cols)
  }
  pkg <- mwperm:::.build_obs_perms_blocks(5, 2L, blocks, ri[keep], ci[keep],
                                          blk, lrow, lcol)
  for (g in pkg) expect_identical(blk[g], blk)        # block-diagonal
  expect_identical(pkg, ref_groups_blocks(ri[keep], ci[keep], blocks, 2L,
                                          5, 1))
})

## ---- incomplete panel and irregular layout (the package's own seams) --------

test_that("incomplete panel: period held fixed inside every block", {
  pn <- make_panel(12, 12, T = 3, beta = 0.1, seed = 89)
  drop <- with_seed(90, sample(unique(paste(pn$i, pn$j)), 10))
  pm <- pn[!(paste(pn$i, pn$j) %in% drop & pn$t > 1), ]
  fit <- mwperm_panel_missing(pm$y, pm$d, x = pm$x, row = pm$i, col = pm$j,
                              time = pm$t, min_block = 3, n_reps = 3,
                              seed = 2, conf_int = FALSE)
  pd <- mwperm:::.panel_missing_design(pm$i, pm$j, pm$t, min_block = 3)
  kd <- pm[pd$idx, ]
  expect_true(all(table(paste(kd$i, kd$j)) == 3))   # every period retained
  groups <- function(r) ref_groups_blocks(dense_id(pm$i)[pd$idx],
                                          dense_id(pm$j)[pd$idx], pd$blocks,
                                          fit$K, 2, r,
                                          slot = dense_id(pm$t)[pd$idx])
  expect_reps_match(fit, kd$y, kd$d, ref_X(kd$x, nrow(kd), time = kd$t),
                    groups)
})

test_that("irregular (Section 6.4): each repetition cuts cells to L0 at random", {
  ly <- make_layout(8, 8, sizes = 3:6, seed = 91, d_cell_constant = TRUE)
  L0 <- 4
  fit <- mwperm_irregular(ly$y, ly$d, x = ly$x, row = ly$i, col = ly$j,
                          rep = ly$l, L0 = L0, min_block = 3, n_reps = 3,
                          seed = 5, conf_int = FALSE)
  pd <- mwperm:::.irregular_design(ly$i, ly$j, ly$l, L0, min_block = 3)
  sub <- ly[pd$idx, ]
  for (r in 1:3) {
    ops <- pd$perm_builder(rep_seed(5, r), fit$K)
    rows <- attr(ops, "rows")
    kd <- sub[rows, ]
    expect_true(all(table(paste(kd$i, kd$j)) == L0))  # exactly L0 per cell
    ref <- ref_procedure1(kd$y, kd$d, ref_X(kd$x, nrow(kd)), unclass(ops))
    expect_gt(ref$gap, 1e-6)
    expect_identical(fit$pvalues_rep[r], ref$pvalue)
  }
  expect_identical(fit$n_obs, as.integer(fit$cells_used * L0))
})
