## find_bicliques(): the fully observed blocks of Definition 1 and Algorithm 2
## (Appendix A) of Guo, Toulis & Wang (2026), which Procedure 2 conditions on.
## Theorem 4 needs, and only needs, blocks that are fully observed and
## disjoint in rows and in columns; the size of the blocks is a power matter
## (Remark 2: an approximate maximum is fine). So validity properties are
## asserted for both methods on random masks, and optimality only where it is
## claimed -- the exact branch-and-bound against brute force on small masks.

## A random mask as a two-column (row, col) matrix of observed cells.
rand_mask <- function(nr, nc, p, seed) {
  A <- with_seed(seed, matrix(stats::runif(nr * nc) < p, nr, nc))
  which(A, arr.ind = TRUE)
}

check_blocks <- function(blocks, cells, min_block) {
  observed <- paste(cells[, 1], cells[, 2])
  used_r <- unlist(lapply(blocks, `[[`, "rows"))
  used_c <- unlist(lapply(blocks, `[[`, "cols"))
  expect_false(anyDuplicated(used_r) > 0)                 # disjoint rows
  expect_false(anyDuplicated(used_c) > 0)                 # disjoint cols
  for (b in blocks) {
    all_cells <- as.vector(outer(b$rows, b$cols, paste))
    expect_true(all(all_cells %in% observed))             # fully observed
    expect_gte(length(b$rows), min_block[1])
    expect_gte(length(b$cols), min_block[2])
  }
}

test_that("blocks are fully observed, disjoint, and respect min_block", {
  for (s in 1:25) {
    nr <- 5 + s %% 7
    nc <- 4 + (3 * s) %% 9
    cells <- rand_mask(nr, nc, p = 0.55 + 0.4 * ((s %% 5) / 5), seed = s)
    for (mb in list(2L, 3L, c(3L, 1L), c(1L, 2L))) {
      mb2 <- if (length(mb) == 1L) rep(max(2L, mb), 2L) else pmax(1L, mb)
      for (method in c("greedy", "exact")) {
        blocks <- suppressWarnings(find_bicliques(
          cells[, 1], cells[, 2], min_block = mb, method = method,
          node_budget = 5e4))
        check_blocks(blocks, cells, mb2)
      }
    }
  }
})

## Largest area over all bicliques, by enumerating every row subset: a row
## set R keeps the columns observed by every row in R.
brute_max_area <- function(A) {
  best <- 0
  nr <- nrow(A)
  for (m in seq_len(2^nr - 1)) {
    R <- which(bitwAnd(m, 2^(seq_len(nr) - 1)) > 0)
    cols <- which(colSums(A[R, , drop = FALSE]) == length(R))
    best <- max(best, length(R) * length(cols))
  }
  best
}

test_that("the exact search finds a maximum-area block (brute force)", {
  for (s in 1:30) {
    nr <- 3 + s %% 5                       # up to 7 rows: 127 subsets
    nc <- 3 + (2 * s) %% 6
    A <- with_seed(100 + s, matrix(stats::runif(nr * nc) < 0.65, nr, nc))
    if (!any(A)) next
    ex <- mwperm:::.max_biclique_exact(A)
    expect_true(ex$exact)
    expect_equal(ex$area, brute_max_area(A))
    expect_true(all(A[ex$rows, ex$cols]))
    ## the greedy heuristic is valid but can be smaller, never larger
    gr <- mwperm:::.grow_biclique(A)
    expect_true(all(A[gr$rows, gr$cols]))
    expect_lte(length(gr$rows) * length(gr$cols), ex$area)
  }
})

test_that("find_bicliques(method = 'exact')'s first block is that maximum", {
  for (s in 1:15) {
    A <- with_seed(200 + s, matrix(stats::runif(36) < 0.75, 6, 6))
    cells <- which(A, arr.ind = TRUE)
    ex <- mwperm:::.max_biclique_exact(A)
    if (min(length(ex$rows), length(ex$cols)) < 2) next
    b1 <- find_bicliques(cells[, 1], cells[, 2], min_block = 2,
                         method = "exact")[[1]]
    expect_equal(length(b1$rows) * length(b1$cols), ex$area)
  }
})

test_that("draft Sec. 6.6: the diagonal-deleted 40 x 40 array gives two 20 x 20 blocks", {
  data(trade_dyadic, package = "mwperm", envir = environment())
  nd <- subset(trade_dyadic, importer != exporter)
  bl <- find_bicliques(as.integer(nd$importer), as.integer(nd$exporter),
                       min_block = 3)
  sz <- sapply(bl, function(b) c(rows = length(b$rows),
                                 cols = length(b$cols)))
  expect_identical(unname(sz), matrix(20L, 2, 2))
  cells <- cbind(as.integer(nd$importer), as.integer(nd$exporter))
  check_blocks(bl, cells, c(3L, 3L))
  ## two blocks of 400 cells: 800 of the 1560 observed cells
  expect_identical(sum(vapply(bl, function(b)
    length(b$rows) * length(b$cols), 0)), 800)
})

test_that("labels come back in the caller's coding, RNG untouched", {
  cells <- rand_mask(8, 8, 0.8, seed = 3)
  lab_r <- c("AUS", "BRA", "CAN", "DEU", "ESP", "FRA", "GBR", "IND")
  set.seed(11)
  before <- .Random.seed
  b_int <- find_bicliques(cells[, 1], cells[, 2], min_block = 2)
  b_chr <- find_bicliques(lab_r[cells[, 1]], cells[, 2], min_block = 2)
  expect_identical(.Random.seed, before)                  # deterministic
  expect_identical(lapply(b_chr, `[[`, "rows"),
                   lapply(b_int, function(b) lab_r[b$rows]))
})

test_that("the exact search falls back to greedy, with a warning, on budget", {
  A <- matrix(TRUE, 30, 30)
  diag(A) <- FALSE
  cells <- which(A, arr.ind = TRUE)
  expect_warning(bl <- find_bicliques(cells[, 1], cells[, 2], min_block = 3,
                                      method = "exact", node_budget = 50),
                 "node budget")
  check_blocks(bl, cells, c(3L, 3L))
})

## ---- a block below the floor does not end the search --------------------------
## Both searches maximise AREA, and a long thin strip has more area than a
## large square block. Until 0.4.4 the peel stopped at the first block below
## min_block, so the strip hid the block: every later block was lost, and
## mwperm_missing() could stop with "no fully observed block" on data that
## had one. Algorithm 2 itself never stops on a small block.
mask_cells <- function(M) {
  idx <- which(M, arr.ind = TRUE)
  list(row = idx[, 1], col = idx[, 2])
}
strip_mask <- function() {              # 2 x 300 strip + a disjoint 20 x 20
  M <- matrix(FALSE, 22, 320)
  M[1:2, 21:320] <- TRUE
  M[3:22, 1:20] <- TRUE
  mask_cells(M)
}
private_rows_mask <- function() {       # 40 rows of 12 private columns + a 3 x 3
  M <- matrix(FALSE, 43, 40 * 12 + 3)
  for (r in 1:40) M[r, (r - 1) * 12 + 1:12] <- TRUE
  M[41:43, 481:483] <- TRUE
  mask_cells(M)
}
sizes_of <- function(bl) vapply(bl, function(b)
  paste0(length(b$rows), "x", length(b$cols)), "")

test_that("a thin strip does not hide a conforming block (both methods)", {
  m <- strip_mask()
  for (meth in c("greedy", "exact")) {
    bl <- find_bicliques(m$row, m$col, min_block = 3L, method = meth)
    expect_identical(sizes_of(bl), "20x20")
    expect_identical(sort(bl[[1]]$rows), as.numeric(3:22))
    expect_identical(sort(bl[[1]]$cols), as.numeric(1:20))
  }
})

test_that("many small private rows do not hide a conforming block", {
  m <- private_rows_mask()
  for (meth in c("greedy", "exact"))
    expect_identical(sizes_of(find_bicliques(m$row, m$col, min_block = 3L,
                                             method = meth)), "3x3")
  N <- length(m$row)
  fit <- with_seed(5, mwperm_missing(stats::rnorm(N), stats::rnorm(N),
                                     row = m$row, col = m$col, min_block = 3L,
                                     seed = 1, conf_int = FALSE))
  expect_identical(fit$n_blocks, 1L)
  expect_identical(fit$n_obs, 9L)
})

test_that("a tie in area between a conforming and a thin block keeps the conforming one", {
  M <- matrix(FALSE, 7, 7)
  M[1:2, 1:2] <- TRUE          # 2 x 2, area 4
  M[3, 3:6] <- TRUE            # 1 x 4, area 4
  M[4:7, 7] <- TRUE            # 4 x 1, area 4
  m <- mask_cells(M)
  for (meth in c("greedy", "exact"))
    expect_true("2x2" %in% sizes_of(find_bicliques(m$row, m$col,
                                                   min_block = 2L,
                                                   method = meth)))
})

test_that("the floor survives a node-budget fallback, and retry_peels = 0 keeps the old stop", {
  m <- strip_mask()
  expect_warning(
    bl <- find_bicliques(m$row, m$col, min_block = 3L, method = "exact",
                         node_budget = 10L),
    "node budget")
  expect_identical(sizes_of(bl), "20x20")
  p <- private_rows_mask()
  expect_length(find_bicliques(p$row, p$col, min_block = 3L,
                               method = "greedy", retry_peels = 0L), 0L)
})

test_that("the floor-constrained exact search is the plain search at the default floor", {
  A <- with_seed(9, matrix(stats::runif(8 * 9) < 0.7, 8, 9))
  expect_identical(mwperm:::.max_biclique_exact(A),
                   mwperm:::.max_biclique_exact(A, mb_r = 1L, mb_c = 1L))
  ## and under a floor it returns a conforming block, or none
  b <- mwperm:::.max_biclique_exact(A, mb_r = 3L, mb_c = 3L)
  expect_true(!length(b$rows) || (length(b$rows) >= 3L && length(b$cols) >= 3L))
  expect_true(all(A[b$rows, b$cols]))
})
