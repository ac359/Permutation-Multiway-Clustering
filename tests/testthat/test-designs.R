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

test_that("layout and irregular: a time-like `rep` draws a warning (0.4.3)", {
  ## Section 6.3 permutes `rep` and Section 6.4 subsamples it at random;
  ## both are exact only for exchangeable replicates, and the authors
  ## (2026-09-30) asked for a warning when the index looks like time.
  ly <- make_layout(4, 4, sizes = 6:8, seed = 85)
  fit <- function(r) mwperm_layout(ly$y, ly$d, row = ly$i, col = ly$j,
                                   rep = r, n_reps = 1, seed = 1,
                                   conf_int = FALSE)
  expect_no_warning(fit(ly$l))                        # 1..n in every cell
  expect_warning(fit(ly$l + 1990), "calendar years")
  expect_warning(fit(ly$l + 3 * (ly$i > 2)), "cells start at different")
  expect_warning(fit(as.Date("2020-01-01") + ly$l), "it holds dates")
  year <- ly$l
  expect_warning(mwperm_layout(ly$y, ly$d, row = ly$i, col = ly$j,
                               rep = year, n_reps = 1, seed = 1,
                               conf_int = FALSE),
                 "its name 'year' is a time word")
  ## the warning changes no number
  expect_identical(suppressWarnings(fit(ly$l + 1990))$pvalues_rep,
                   fit(ly$l)$pvalues_rep)
  li <- make_layout(8, 8, sizes = 3:6, seed = 107, d_cell_constant = TRUE)
  expect_warning(
    mwperm_irregular(li$y, li$d, row = li$i, col = li$j, rep = li$l + 2000,
                     L0 = 3, min_block = 3, n_reps = 2, seed = 1,
                     conf_int = FALSE),
    "mwperm_panel_missing(time = <rep>, L0 = 3)", fixed = TRUE)
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
  ## a K no block can carry is refused
  expect_error(mwperm_missing(md$y, md$d, x = md$x, row = md$i, col = md$j,
                              min_block = 3, K = max(sides), seed = 1),
               "largest fully observed block")
})

test_that("missing: blocks smaller than K + 1 are kept and held fixed (0.4.3)", {
  ## The method's authors (2026-09-30): it is valid to keep blocks smaller
  ## than K + 1 and leave them unpermuted. An explicit K above the default
  ## does that; the reference applies the identity to those blocks.
  set.seed(90)
  grp <- c(rep(1L, 12), rep(2L, 5), rep(3L, 3))
  g <- expand.grid(i = seq_len(20), j = seq_len(20))
  g <- g[grp[g$i] == grp[g$j], ]
  g$d <- stats::rnorm(nrow(g)); g$x <- stats::rnorm(nrow(g))
  g$y <- 0.2 * g$d + stats::rnorm(20)[g$i] + stats::rnorm(nrow(g))
  fit <- mwperm_missing(g$y, g$d, x = g$x, row = g$i, col = g$j,
                        min_block = 3, K = 9, n_reps = 3, seed = 7,
                        conf_int = FALSE)
  expect_identical(fit$K, 9L)
  expect_identical(fit$cells_used, nrow(g))           # nothing discarded
  expect_true(any(grepl("held fixed", fit$note, fixed = TRUE)))
  blocks <- find_bicliques(dense_id(g$i), dense_id(g$j), min_block = 3)
  groups <- function(r) ref_groups_blocks(dense_id(g$i), dense_id(g$j),
                                          blocks, 9L, 7, r)
  expect_reps_match(fit, g$y, g$d, ref_X(g$x, nrow(g)), groups)
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

## ---- layout: K above the smallest cell (Section 6.3 with Algorithm 1) --------
## One K for all cells; a cell with fewer than K + 1 replicates is left
## unpermuted and stays in the regression. Before 0.4.5 the smallest cell
## capped K, so one cell of two replicates made the smallest p-value 1/2 and a
## singleton cell was an error.
held_layout <- function(sizes = c(2L, 25L, 30L, 22L, 1L, 40L, 20L, 21L, 26L),
                        seed = 31) {
  with_seed(seed, {
    cells <- expand.grid(i = 1:3, j = 1:3)
    g <- data.frame(i = rep(cells$i, sizes), j = rep(cells$j, sizes),
                    l = unlist(lapply(sizes, seq_len)))
    cid <- rep(seq_len(nrow(cells)), sizes)
    g$d <- stats::rnorm(nrow(g))
    g$x <- stats::rnorm(nrow(g))
    g$y <- 0.5 * g$x + stats::rnorm(nrow(cells))[cid] + stats::rnorm(nrow(g))
    g
  })
}

test_that("layout: an explicit K above the smallest cell holds small cells fixed", {
  ly <- held_layout()
  fit <- mwperm_layout(ly$y, ly$d, x = ly$x, row = ly$i, col = ly$j,
                       rep = ly$l, K = 19, seed = 1, conf_int = FALSE)
  expect_identical(fit$K, 19L)
  expect_equal(fit$p_floor, 0.05)
  expect_identical(fit$n_obs, nrow(ly))                  # nothing is dropped
  expect_true(any(grepl("2 cells with fewer than K + 1 = 20 replicates",
                        fit$note, fixed = TRUE)))
  ## the groups: the identity for the cells of 2 and 1, Algorithm 1 elsewhere
  cell <- dense_id(interaction(dense_id(ly$i), dense_id(ly$j), drop = TRUE))
  ell <- tabulate(cell, max(cell))
  G <- mwperm:::.layout_groups(ell, 19L, rep_seed(1, 1), 1000)
  for (c in seq_along(ell)) {
    if (ell[c] < 20L)
      expect_true(all(vapply(G[[c]], identical, NA, seq_len(ell[c]))))
    else
      expect_identical(G[[c]], build_perm_set(ell[c], 19L,
                                              seed = sub_seed(rep_seed(1, 1), c)))
  }
  ## the p-value is the naive Procedure 1 on the literal construction
  ref <- ref_procedure1(ly$y, ly$d, ref_X(ly$x, nrow(ly)),
                        ref_groups_layout(ly$i, ly$j, ly$l, 19L, 1, 1))
  expect_equal(fit$pvalue, ref$pvalue, tolerance = 1e-12)
  ## seeded: the same call gives the same object
  fit2 <- mwperm_layout(ly$y, ly$d, x = ly$x, row = ly$i, col = ly$j,
                        rep = ly$l, K = 19, seed = 1, conf_int = FALSE)
  expect_identical(fit$pvalues_rep, fit2$pvalues_rep)
})

test_that("layout: the default K permutes every cell with two or more replicates", {
  ly <- held_layout()
  expect_no_error(
    fit <- mwperm_layout(ly$y, ly$d, x = ly$x, row = ly$i, col = ly$j,
                         rep = ly$l, seed = 1, conf_int = FALSE))
  expect_identical(fit$K, 1L)             # the smallest cell with >= 2 replicates
  expect_true(any(grepl("1 cell with fewer than K + 1 = 2 replicates",
                        fit$note, fixed = TRUE)))
  expect_true(any(grepl("A larger `K` is available", fit$note, fixed = TRUE)))
  ## without a singleton nothing changes: the default is the smallest cell - 1
  l2 <- make_layout(4, 4, sizes = 6:8, seed = 57)
  f2 <- mwperm_layout(l2$y, l2$d, x = l2$x, row = l2$i, col = l2$j,
                      rep = l2$l, seed = 2, conf_int = FALSE)
  expect_identical(f2$K, 5L)
  expect_false(any(grepl("with fewer than K + 1 =", f2$note, fixed = TRUE)))
})

test_that("layout: a K no cell can carry, and a layout with no replication, stop", {
  ly <- held_layout()
  expect_error(
    mwperm_layout(ly$y, ly$d, row = ly$i, col = ly$j, rep = ly$l, K = 60),
    "largest cell holds only 40 replicates")
  one <- expand.grid(i = 1:4, j = 1:4)
  expect_error(
    mwperm_layout(seq_len(16) + 0, rev(seq_len(16)) + 0, row = one$i,
                  col = one$j),
    "No cell holds two or more observations")
})

test_that("three-way: one small dimension draws the note that names the panel route", {
  ## Algorithm 1 with a common K leaves a dimension with fewer than K + 1
  ## levels fixed; that test is the panel design with that dimension as time,
  ## valid under three-way exchangeability. A 25 x 22 x 6 array caps K at 5.
  g <- with_seed(41, {
    g <- expand.grid(a = 1:25, b = 1:22, c = 1:6)
    g$d <- stats::rnorm(nrow(g)); g$y <- stats::rnorm(nrow(g)); g
  })
  f <- mwperm_threeway(g$y, g$d, id1 = g$a, id2 = g$b, id3 = g$c, seed = 3,
                       conf_int = FALSE)
  expect_identical(f$K, 5L)
  expect_true(any(grepl("Only `id3` is too small", f$note, fixed = TRUE)))
  expect_true(any(grepl("mwperm_panel() with `time = ` that index reaches K = 21",
                        f$note, fixed = TRUE)))
  ## and that route is the literal construction: Algorithm 1 three times with
  ## K = 21, the third dimension (6 levels) left fixed
  p <- mwperm_panel(g$y, g$d, row = g$a, col = g$b, time = g$c,
                    time_fe = FALSE, K = 21, seed = 3, conf_int = FALSE)
  expect_identical(p$K, 21L)
  ## no such note when every dimension is large enough, or when two are small
  g2 <- with_seed(42, {
    g <- expand.grid(a = 1:6, b = 1:7, c = 1:5)
    g$d <- stats::rnorm(nrow(g)); g$y <- stats::rnorm(nrow(g)); g
  })
  f2 <- mwperm_threeway(g2$y, g2$d, id1 = g2$a, id2 = g2$b, id3 = g2$c,
                        seed = 3, conf_int = FALSE)
  expect_false(any(grepl("is too small here", f2$note, fixed = TRUE)))
})

test_that("a coarse fit with several repetitions does not call the median exact", {
  dy <- make_dyadic(6, 6, seed = 71)
  f1 <- mwperm_dyadic(dy$y, dy$d, x = dy$x, row = dy$i, col = dy$j, seed = 1)
  expect_true(any(grepl("remains exact", f1$note, fixed = TRUE)))   # one run
  f5 <- mwperm_dyadic(dy$y, dy$d, x = dy$x, row = dy$i, col = dy$j, seed = 1,
                      n_reps = 5)
  expect_false(any(grepl("reported above is unaffected and remains exact",
                         f5$note, fixed = TRUE)))
  expect_true(any(grepl("their median over 5 repetitions, which the theorem does not cover",
                        f5$note, fixed = TRUE)))
  f5b <- mwperm_dyadic(dy$y, dy$d, x = dy$x, row = dy$i, col = dy$j, seed = 1,
                       n_reps = 5, aggregate = "median2")
  expect_true(any(grepl("twice their median is valid", f5b$note, fixed = TRUE)))
})
