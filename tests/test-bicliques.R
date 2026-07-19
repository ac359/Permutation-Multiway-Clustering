## find_bicliques() contract tests (audit Phase 5, brief section 5 item 7 /
## CLAUDE.md section 5.6). Base-R stopifnot style; fast.
## The contract mwperm_missing() relies on: every returned block is FULLY
## OBSERVED; blocks are DISJOINT in both rows and columns; min_block honoured
## (floored at 2); deterministic (no RNG, caller's stream untouched); the
## exact method's node-budget fallback warns and still returns valid blocks.
library(mwperm)

msg_of <- function(expr)
  tryCatch({
    expr
    NA_character_
  }, error = function(e) conditionMessage(e))

## contract checker: TRUE iff all invariants hold for one block list
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

## ---- 1. random-mask sweep, both methods -------------------------------------
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

## ---- 2. determinism + RNG hygiene -------------------------------------------
g <- expand.grid(i = 1:10, j = 1:10)
set.seed(5)
g <- g[runif(nrow(g)) < .7, ]
b1 <- find_bicliques(g$i, g$j, min_block = 3)
set.seed(999)
b2 <- find_bicliques(g$i, g$j, min_block = 3)
stopifnot(identical(b1, b2))
invisible(runif(1))                       # ensure .Random.seed exists
rs <- .Random.seed
invisible(find_bicliques(g$i, g$j, min_block = 3))
stopifnot(identical(rs, .Random.seed))

## ---- 3. node-budget fallback: warning + valid blocks ------------------------
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

## ---- 4. min_block floor + label restoration ---------------------------------
gc2 <- expand.grid(i = c("IT", "FR", "DE", "US"), j = c("w", "x", "y", "z"),
                   stringsAsFactors = FALSE)
b_floor <- find_bicliques(gc2$i, gc2$j, min_block = -5)   # floored to 2
stopifnot(length(b_floor) >= 1L,
          min(vapply(b_floor, function(b) min(length(b$rows), length(b$cols)),
                     integer(1))) >= 2L,
          all(unlist(lapply(b_floor, `[[`, "rows")) %in% gc2$i))
bln <- find_bicliques(c(10, 20, 35), c(7, 7, 7), min_block = 2)  # 3x1: no block
stopifnot(length(bln) == 0L)
gnum <- expand.grid(i = c(10, 20, 35), j = c(7, 9, 11))
bln2 <- find_bicliques(gnum$i, gnum$j, min_block = 2)
stopifnot(is.numeric(bln2[[1]]$rows), all(bln2[[1]]$rows %in% c(10, 20, 35)))

## ---- 5. degenerate masks + the mwperm_missing() no-block error --------------
stopifnot(length(find_bicliques(1:5, 1:5,
                                min_block = 2)) == 0L)  # diagonal only
m <- msg_of(mwperm_missing(rnorm(5), rnorm(5), row = 1:5, col = 1:5, seed = 1))
stopifnot(!is.na(m), grepl("min_block", m))

## ---- 6. mwperm_missing() uses exactly the blocks find_bicliques() returns ---
set.seed(7)
g8 <- expand.grid(i = 1:8, j = 1:8)
g8 <- g8[g8$i != g8$j, ]
y8 <- rnorm(nrow(g8))
d8 <- rnorm(nrow(g8))
bl8 <- find_bicliques(g8$i, g8$j, min_block = 3)
fit <- mwperm_missing(y8, d8, row = g8$i, col = g8$j, min_block = 3, seed = 1,
                      conf_int = FALSE)
in_bl <- rep(FALSE, nrow(g8))
for (b in bl8) in_bl <- in_bl | (g8$i %in% b$rows & g8$j %in% b$cols)
stopifnot(fit$n_blocks == length(bl8), fit$cells_used == sum(in_bl),
          fit$K == min(vapply(bl8, function(b)
            min(length(b$rows), length(b$cols)), integer(1))) - 1L)

cat("test-bicliques.R: all assertions passed\n")
