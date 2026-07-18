## One-sided permutation for the missing-cells design (permute = "rows"/
## "cols"): the builder must apply a genuine subgroup -- the held-fixed
## dimension untouched, the permuted one a closed cyclic group -- and the
## front end must unlock K via the permuted side only. Base-R stopifnot.
library(mwperm)

msg_of <- function(expr) tryCatch({ expr; NA_character_ },
                                  error = function(e) conditionMessage(e))
warns_of <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
  })
  w
}

## ---- 1. builder contract: rows-only is a subgroup, columns untouched ---------
blocks <- list(list(rows = 1:6, cols = 1:4), list(rows = 7:10, cols = 5:6))
cells <- do.call(rbind, lapply(seq_along(blocks), function(q) {
  b <- blocks[[q]]
  g <- expand.grid(r = b$rows, c = b$cols)
  data.frame(ri = g$r, ci = g$c, blk = q,
             lrow = match(g$r, b$rows), lcol = match(g$c, b$cols))
}))
K <- 3L
ops <- mwperm:::.build_obs_perms_blocks(11L, K, blocks,
                                        ri = cells$ri, ci = cells$ci,
                                        blk = cells$blk, lrow = cells$lrow,
                                        lcol = cells$lcol, permute = "rows")
N <- nrow(cells)
stopifnot(length(ops) == K + 1L, identical(ops[[1L]], seq_len(N)))
for (g in ops) stopifnot(identical(sort(g), seq_len(N)))     # permutations
for (g in ops) stopifnot(identical(cells$ci[g], cells$ci))   # cols fixed
stopifnot(any(cells$ri[ops[[2L]]] != cells$ri))              # rows really move
## closure under composition (cyclic group)
keys <- vapply(ops, paste, "", collapse = ",")
for (a in seq_along(ops)) for (b in seq_along(ops))
  stopifnot(paste(ops[[a]][ops[[b]]], collapse = ",") %in% keys)
## cols-only mirrors: rows fixed, columns move
opc <- mwperm:::.build_obs_perms_blocks(11L, 1L, blocks,
                                        ri = cells$ri, ci = cells$ci,
                                        blk = cells$blk, lrow = cells$lrow,
                                        lcol = cells$lcol, permute = "cols")
for (g in opc) stopifnot(identical(cells$ri[g], cells$ri))
stopifnot(any(cells$ci[opc[[2L]]] != cells$ci))

## ---- 2. find_bicliques: length-2 min_block contract ---------------------------
## a 4x1 strip: usable for one-sided permutation, invisible to the scalar floor
strip <- data.frame(i = 1:4, j = rep(1L, 4))
b31 <- find_bicliques(strip$i, strip$j, min_block = c(3, 1))
stopifnot(length(b31) == 1L, length(b31[[1L]]$rows) == 4L,
          length(b31[[1L]]$cols) == 1L)
stopifnot(length(find_bicliques(strip$i, strip$j, min_block = 2)) == 0L)
## scalar behaviour is unchanged: min_block = m == c(m, m)
g2 <- expand.grid(i = 1:6, j = 1:6); gm <- g2[g2$i != g2$j, ]
stopifnot(identical(find_bicliques(gm$i, gm$j, min_block = 2),
                    find_bicliques(gm$i, gm$j, min_block = c(2, 2))))
stopifnot(grepl("at least one side",
                msg_of(find_bicliques(strip$i, strip$j, min_block = c(1, 1)))),
          grepl("one integer",
                msg_of(find_bicliques(strip$i, strip$j, min_block = c(NA, 2)))))
## asymmetric floor vs area maximisation: the 6x7 dense block dominates area,
## but only the 10x1 strip satisfies c(8, 1) -- the constrained retry finds it
wide <- rbind(expand.grid(i = 1:6, j = 2:7), data.frame(i = 1:10, j = 1))
b81 <- find_bicliques(wide$i, wide$j, min_block = c(8, 1))
stopifnot(length(b81) == 1L, identical(b81[[1L]]$rows, as.numeric(1:10)),
          identical(b81[[1L]]$cols, 1))
## transposed: c(1, 8) on the transposed mask finds the 1x10 wide block
b18 <- find_bicliques(wide$j, wide$i, min_block = c(1, 8))
stopifnot(length(b18) == 1L, length(b18[[1L]]$cols) == 10L)

## ---- 3. K unlock: the permuted side alone caps the group order ----------------
## single fully observed column, 30 rows: two-sided has no usable block at
## all; rows-only reaches K = 29 (resolution 1/30 <= 0.05)
set.seed(7)
n1 <- 30L
eta <- rnorm(n1)                                   # row effects (exchangeable)
d1 <- rnorm(n1)
y1 <- 2 * d1 + eta + rnorm(n1, sd = 0.3)
f1 <- mwperm_missing(y1, d1, row = seq_len(n1), col = rep(1L, n1),
                     permute = "rows", seed = 3, n_reps = 2)
stopifnot(f1$K == 29L, f1$type == "missing (bicliques, rows-only)",
          f1$cells_used == n1, f1$pvalue <= 0.05,   # strong true effect
          !is.null(f1$conf_int), all(is.finite(f1$conf_int)))
stopifnot(grepl("No fully observed block",
                msg_of(mwperm_missing(y1, d1, row = seq_len(n1),
                                      col = rep(1L, n1), seed = 3))))
## seeded determinism
f1b <- mwperm_missing(y1, d1, row = seq_len(n1), col = rep(1L, n1),
                      permute = "rows", seed = 3, n_reps = 2)
stopifnot(identical(f1$pvalue, f1b$pvalue),
          identical(f1$estimate, f1b$estimate),
          identical(f1$conf_int, f1b$conf_int))

## ---- 4. resolution note names the permuted side -------------------------------
## rows-only on short blocks: the cap message must exist and cite the side
gs <- expand.grid(i = 1:3, j = 1:8)                # 3 rows x 8 cols, complete
fs <- mwperm_missing(rnorm(24), rnorm(24), row = gs$i, col = gs$j,
                     permute = "rows", min_block = 2, seed = 1, n_reps = 1,
                     conf_int = FALSE)
stopifnot(fs$K == 2L, any(grepl("permuted side", fs$note)))
## same data, cols-only: the 8-column side unlocks K = 7
fc <- mwperm_missing(rnorm(24), rnorm(24), row = gs$i, col = gs$j,
                     permute = "cols", min_block = 2, seed = 1, n_reps = 1,
                     conf_int = FALSE)
stopifnot(fc$K == 7L, fc$type == "missing (bicliques, cols-only)")

## ---- 5. mwperm() passthrough + wrong-design warning ---------------------------
ym <- rnorm(nrow(gm)); dm <- rnorm(nrow(gm))
dir5 <- mwperm_missing(ym, dm, row = gm$i, col = gm$j, permute = "rows",
                       min_block = 3, seed = 5, n_reps = 2, conf_int = FALSE)
via5 <- mwperm(y = ym, d = dm, index = list(i = gm$i, j = gm$j),
               permute = "rows", min_block = 3, seed = 5, n_reps = 2,
               conf_int = FALSE, verbose = FALSE)
stopifnot(identical(via5$auto$design, "missing"),
          identical(dir5$pvalue, via5$pvalue),
          identical(dir5$estimate, via5$estimate),
          identical(dir5$K, via5$K), identical(dir5$type, via5$type))
w5 <- warns_of(mwperm(y = rnorm(36), d = rnorm(36),
                      index = list(i = g2$i, j = g2$j), permute = "rows",
                      seed = 1, n_reps = 1, conf_int = FALSE,
                      verbose = FALSE))
stopifnot(any(grepl("`permute` applies to the missing design only", w5,
                    fixed = TRUE)))

cat("test-permute.R: all assertions passed\n")
