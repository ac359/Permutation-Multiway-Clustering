## mwperm_missing() and find_bicliques() -- incomplete arrays, Procedure 2.
##
## When cells are missing there is no relabelling of the whole array that maps
## observed cells to observed cells, so the test restricts to fully observed
## rectangular blocks (bicliques) that are DISJOINT in both rows and columns,
## permutes inside each block with a common group order, and pools the residual
## statistics (Guo, Toulis & Wang 2026, Section 5). Discarded cells are the
## price of exact validity; a sub-maximal block costs power, never validity.
##
## find_bicliques() is exported and tested here alongside the front end that
## consumes it, because the front end's correctness rests on its contract:
## every block fully observed, blocks disjoint in both dimensions, no RNG.
##
## Cross-cutting contracts live elsewhere; see tests/README.md. The block
## gather-vector algebra is in `lower-level-tests/test-obsperms.R`.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

## TRUE iff a block list satisfies every invariant mwperm_missing() relies on
blocks_ok <- function(blocks, ri, ci, min_block) {
  obs <- paste(ri, ci)
  seen_r <- c()
  seen_c <- c()
  for (b in blocks) {
    if (!all(as.vector(outer(b$rows, b$cols, paste)) %in% obs)) return(FALSE)
    if (min(length(b$rows), length(b$cols)) < max(2L, min_block)) return(FALSE)
    if (length(intersect(b$rows, seen_r)) || length(intersect(b$cols, seen_c)))
      return(FALSE)
    seen_r <- c(seen_r, b$rows)
    seen_c <- c(seen_c, b$cols)
  }
  TRUE
}

## ---- 1. find_bicliques(): the contract, over random masks -----------------
set.seed(42)
for (m in 1:15) {
  nr <- sample(5:10, 1)
  nc <- sample(5:10, 1)
  g <- expand.grid(i = seq_len(nr), j = seq_len(nc))
  g <- g[runif(nrow(g)) < runif(1, 0.4, 0.95), ]
  if (nrow(g) < 4L) next
  mb <- sample(2:3, 1)
  for (meth in c("greedy", "exact")) {
    bl <- suppressWarnings(find_bicliques(g$i, g$j, min_block = mb,
                                          method = meth))
    stopifnot(blocks_ok(bl, g$i, g$j, mb))
  }
}

## ---- 2. find_bicliques() is deterministic and RNG-free --------------------
## A seeded fit must be reproducible, so the block search may not consume the
## caller's random stream -- nor depend on it.
g <- expand.grid(i = 1:10, j = 1:10)
set.seed(5)
g <- g[runif(nrow(g)) < .7, ]
b1 <- find_bicliques(g$i, g$j, min_block = 3)
set.seed(999)
b2 <- find_bicliques(g$i, g$j, min_block = 3)
stopifnot(identical(b1, b2))
invisible(runif(1))                                # ensure .Random.seed exists
rs <- .Random.seed
invisible(find_bicliques(g$i, g$j, min_block = 3))
stopifnot(identical(rs, .Random.seed))

## ---- 3. the exact search's node-budget fallback ---------------------------
## Maximum biclique is NP-hard. On hitting the budget the branch and bound
## falls back to the greedy block: still valid, just possibly sub-maximal, and
## it must say so rather than fail or return silently worse blocks.
g20 <- expand.grid(i = 1:20, j = 1:20)
g20 <- g20[g20$i != g20$j, ]
w <- NA_character_
bl20 <- withCallingHandlers(
  find_bicliques(g20$i, g20$j, min_block = 3, method = "exact",
                 node_budget = 50L),
  warning = function(x) {
    w <<- conditionMessage(x)
    invokeRestart("muffleWarning")
  })
stopifnot(!is.na(w), grepl("node budget", w),
          blocks_ok(bl20, g20$i, g20$j, 3L))

## ---- 4. min_block: the floor, asymmetric floors, and labels ---------------
gc2 <- expand.grid(i = c("IT", "FR", "DE", "US"), j = c("w", "x", "y", "z"),
                   stringsAsFactors = FALSE)
b_floor <- find_bicliques(gc2$i, gc2$j, min_block = -5)      # floored to 2
stopifnot(length(b_floor) >= 1L,
          min(vapply(b_floor, function(b) min(length(b$rows), length(b$cols)),
                     integer(1))) >= 2L,
          all(unlist(lapply(b_floor, `[[`, "rows")) %in% gc2$i))
## the original labels come back, in their original type
bln2 <- find_bicliques(expand.grid(i = c(10, 20, 35), j = c(7, 9, 11))$i,
                       expand.grid(i = c(10, 20, 35), j = c(7, 9, 11))$j,
                       min_block = 2)
stopifnot(is.numeric(bln2[[1]]$rows), all(bln2[[1]]$rows %in% c(10, 20, 35)))
## a length-2 min_block sets the two sides separately: a 4 x 1 strip is usable
## for one-sided permutation but invisible to the scalar floor
strip <- data.frame(i = 1:4, j = rep(1L, 4))
b31 <- find_bicliques(strip$i, strip$j, min_block = c(3, 1))
stopifnot(length(b31) == 1L, length(b31[[1L]]$rows) == 4L,
          length(b31[[1L]]$cols) == 1L,
          length(find_bicliques(strip$i, strip$j, min_block = 2)) == 0L)
## scalar behaviour is unchanged: min_block = m is exactly c(m, m)
g6 <- expand.grid(i = 1:6, j = 1:6)
gm6 <- g6[g6$i != g6$j, ]
stopifnot(identical(find_bicliques(gm6$i, gm6$j, min_block = 2),
                    find_bicliques(gm6$i, gm6$j, min_block = c(2, 2))))
expect_err(find_bicliques(strip$i, strip$j, min_block = c(1, 1)),
           "at least one side")
expect_err(find_bicliques(strip$i, strip$j, min_block = c(NA, 2)),
           "one integer")
## asymmetric floor vs area maximisation: the 6 x 7 dense block dominates on
## area, but only the 10 x 1 strip satisfies c(8, 1) -- the constrained retry
## has to find it
wide <- rbind(expand.grid(i = 1:6, j = 2:7), data.frame(i = 1:10, j = 1))
b81 <- find_bicliques(wide$i, wide$j, min_block = c(8, 1))
stopifnot(length(b81) == 1L, identical(b81[[1L]]$rows, as.numeric(1:10)),
          identical(b81[[1L]]$cols, 1))
b18 <- find_bicliques(wide$j, wide$i, min_block = c(1, 8))   # transposed
stopifnot(length(b18) == 1L, length(b18[[1L]]$cols) == 10L)

## ---- 5. degenerate masks, and the no-usable-block error -------------------
stopifnot(length(find_bicliques(1:5, 1:5, min_block = 2)) == 0L)  # diagonal
expect_err(mwperm_missing(rnorm(5), rnorm(5), row = 1:5, col = 1:5, seed = 1),
           "min_block")
stopifnot(grepl("No fully observed block",
                msg_of(mwperm_missing(rnorm(5), rnorm(5), row = 1:5,
                                      col = 1:5, seed = 1))))

## ---- 6. the fit uses exactly the blocks find_bicliques() returns ----------
set.seed(7)
g8 <- expand.grid(i = 1:8, j = 1:8)
g8 <- g8[g8$i != g8$j, ]
N8 <- nrow(g8)
y8 <- rnorm(N8)
d8 <- rnorm(N8)
bl8 <- find_bicliques(g8$i, g8$j, min_block = 3)
fit8 <- mwperm_missing(y8, d8, row = g8$i, col = g8$j, min_block = 3, seed = 1,
                       conf_int = FALSE)
in_bl <- rep(FALSE, N8)
for (b in bl8) in_bl <- in_bl | (g8$i %in% b$rows & g8$j %in% b$cols)
stopifnot(identical(fit8$type, "missing (bicliques)"),
          fit8$n_blocks == length(bl8),
          fit8$cells_used == sum(in_bl),
          fit8$cells_total == N8,
          fit8$K == min(vapply(bl8, function(b)
            min(length(b$rows), length(b$cols)), integer(1))) - 1L)
## the smallest block caps K, and the note names that as the cause -- with an
## attainable alpha there is nothing to explain and the note is absent
stopifnot(1 / (fit8$K + 1) > 0.05,            # the fixture really is capped
          any(grepl("smallest selected block", fit8$note)))
fit8a <- mwperm_missing(y8, d8, row = g8$i, col = g8$j, min_block = 3, seed = 1,
                        conf_int = FALSE, alpha = 0.4)
stopifnot(!any(grepl("smallest selected block", fit8a$note)))

## ---- 7. validation runs on the FULL data, before cells are discarded ------
## The observed cells are the user's data contract: a non-finite value must be
## an error whether or not the biclique step would have thrown that cell away.
i_kept <- which(g8$i %in% bl8[[1]]$rows & g8$j %in% bl8[[1]]$cols)[1L]
y_bad <- y8
y_bad[i_kept] <- NA
expect_err(mwperm_missing(y_bad, d8, row = g8$i, col = g8$j, min_block = 3),
           "`y`")
if (any(!in_bl)) {
  y_bad2 <- y8
  y_bad2[which(!in_bl)[1L]] <- NA
  expect_err(mwperm_missing(y_bad2, d8, row = g8$i, col = g8$j, min_block = 3),
             "`y`")
}
expect_err(mwperm_missing(c(y8, 1), c(d8, 1), row = c(g8$i, 1),
                          col = c(g8$j, 2)),
           "one observation per (row, col) cell")

## ---- 8. the seeded anchor on the shipped data -----------------------------
data(trade_dyadic)
set.seed(1)
dm <- trade_dyadic[trade_dyadic$importer != trade_dyadic$exporter, ]
dm <- dm[sample(nrow(dm), round(0.9 * nrow(dm))), ]
fit <- with(dm, mwperm_missing(y = log_trade, d = log_dist,
                               x = cbind(log_gdp_i, log_gdp_j),
                               row = importer, col = exporter,
                               min_block = 3, seed = 1))
stopifnot(abs(fit$estimate - (-0.804887)) < 1e-6,
          identical(fit$pvalue, 0.2),
          fit$K == 4L, fit$n_blocks == 4L,
          fit$cells_used == 433L, fit$cells_total == 1404L,
          is.null(fit$conf_int))                   # 1/(K+1) = 0.2 > alpha

## ---- 9. one-sided permutation unlocks K on lopsided blocks ----------------
## permute = "rows"/"cols" applies a genuine subgroup: the same exchangeability
## assumption covers it, and K is capped by the PERMUTED side alone. A single
## fully observed column has no two-sided block at all, yet 30 exchangeable
## rows give K = 29 and so a usable 95% interval.
set.seed(7)
n1 <- 30L
d1 <- rnorm(n1)
y1 <- 2 * d1 + rnorm(n1) + rnorm(n1, sd = 0.3)
f_rows <- mwperm_missing(y1, d1, row = seq_len(n1), col = rep(1L, n1),
                         permute = "rows", seed = 3, n_reps = 2)
stopifnot(f_rows$K == 29L,
          identical(f_rows$type, "missing (bicliques, rows-only)"),
          f_rows$cells_used == n1,
          f_rows$pvalue <= 0.05,                   # a strong true effect
          !is.null(f_rows$conf_int), all(is.finite(f_rows$conf_int)))
stopifnot(grepl("No fully observed block",
                msg_of(mwperm_missing(y1, d1, row = seq_len(n1),
                                      col = rep(1L, n1), seed = 3))))
stopifnot(isTRUE(same_fit(
  f_rows, mwperm_missing(y1, d1, row = seq_len(n1), col = rep(1L, n1),
                         permute = "rows", seed = 3, n_reps = 2))))

## the resolution note cites the permuted side, and the other side unlocks a
## different K on the very same data
gs <- expand.grid(i = 1:3, j = 1:8)                # 3 rows x 8 cols, complete
set.seed(9)
ys <- rnorm(24)
ds <- rnorm(24)
f_r <- mwperm_missing(ys, ds, row = gs$i, col = gs$j, permute = "rows",
                      min_block = 2, seed = 1, n_reps = 1, conf_int = FALSE)
f_c <- mwperm_missing(ys, ds, row = gs$i, col = gs$j, permute = "cols",
                      min_block = 2, seed = 1, n_reps = 1, conf_int = FALSE)
stopifnot(f_r$K == 2L, any(grepl("permuted side", f_r$note)),
          f_c$K == 7L,
          identical(f_c$type, "missing (bicliques, cols-only)"))

## mwperm() forwards `permute` unchanged
set.seed(11)
ym <- rnorm(nrow(gm6))
dmm <- rnorm(nrow(gm6))
direct <- mwperm_missing(ym, dmm, row = gm6$i, col = gm6$j, permute = "rows",
                         min_block = 3, seed = 5, n_reps = 2, conf_int = FALSE)
viad <- mwperm(y = ym, d = dmm, index = list(i = gm6$i, j = gm6$j),
               permute = "rows", min_block = 3, seed = 5, n_reps = 2,
               conf_int = FALSE, verbose = FALSE)
stopifnot(identical(viad$auto$design, "missing"),
          identical(direct$pvalue, viad$pvalue),
          identical(direct$estimate, viad$estimate),
          identical(direct$K, viad$K), identical(direct$type, viad$type))

passed("test-missing.R")
