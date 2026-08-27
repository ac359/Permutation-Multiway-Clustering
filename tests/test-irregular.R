## mwperm_irregular(): the Section 6.4 procedure.
##
## The design claim is structural: on a complete balanced array with L0 equal to
## the number of within-cell slots, Section 6.4 degenerates to the panel
## construction of Section 6.2 -- the mask is all ones, the biclique search
## returns the whole array as one block, nothing is deleted, and Procedure 2's
## row/column groups act identically in every slot with the slot held fixed.
## That is exactly what mwperm_panel() does with time. This file pins the
## gather vectors themselves, not just a p-value.
##
## Reaches internals via ::: -- run against a FRESHLY INSTALLED package.
library(mwperm)

ok <- function(...) stopifnot(...)
err <- function(expr, pat) {
  m <- tryCatch({ expr; NA_character_ },
                error = function(e) conditionMessage(e))
  ok(!is.na(m), grepl(pat, m, fixed = TRUE))
}

data(trade_panel)
tp <- trade_panel
ri <- as.integer(factor(tp$importer))
ci <- as.integer(factor(tp$exporter))
ti <- as.integer(factor(tp$year))
m <- max(ri); n <- max(ci); TT <- max(ti)

## ---- 1. same permutations as the panel construction ------------------------
## The two front ends draw their groups from different sub-seed offsets (1, 2
## for the panel; 4q-1, 4q for the block builder), so the groups are matched by
## hand here and the resulting observation gather vectors compared directly.
K <- m - 1L
rs <- 1L
Gr <- build_perm_set(m, K, seed = mwperm:::.sub_seed(rs, 3L))
Gc <- build_perm_set(n, K, seed = mwperm:::.sub_seed(rs, 4L))
A <- mwperm:::.build_obs_perms(cbind(ri, ci, ti), list(Gr, Gc, NULL),
                               design = "panel")
B <- mwperm:::.build_obs_perms_blocks(
  rs, K, list(list(rows = seq_len(m), cols = seq_len(n))),
  ri = ri, ci = ci, blk = rep(1L, length(ri)), lrow = ri, lcol = ci, slot = ti)
ok(identical(A, B))

## ---- 2. same data retained, hence same OLS reference -----------------------
## L0 = number of periods keeps every observation, so the two fits see exactly
## the same design matrix and must agree on estimate, naive SE, N and K.
fp <- with(tp, mwperm_panel(y = log_trade, d = fta,
                            x = cbind(log_gdp_i, log_gdp_j),
                            row = importer, col = exporter, time = year,
                            time_fe = FALSE, conf_int = FALSE, n_reps = 2L,
                            seed = 1))
fi <- with(tp, mwperm_irregular(y = log_trade, d = fta,
                                x = cbind(log_gdp_i, log_gdp_j),
                                row = importer, col = exporter, rep = year,
                                L0 = TT, min_block = 2L, conf_int = FALSE,
                                n_reps = 2L, seed = 1))
ok(identical(unname(fp$estimate), unname(fi$estimate)),
   identical(unname(fp$se_naive), unname(fi$se_naive)),
   identical(fp$n_obs, fi$n_obs), identical(fp$K, fi$K),
   fi$n_blocks == 1L, fi$cells_used == m * n, fi$L0 == TT)

## ---- 3. it works where the within-cell test has none -----------------------
## d constant within every cell: mwperm_layout() warns and returns p = 1 by
## construction; mwperm_irregular() permutes the cells instead and can reject.
set.seed(11)
dat <- do.call(rbind, lapply(1:8, function(i)
  do.call(rbind, lapply(1:8, function(j) {
    L <- sample(c(0L, 2L, 4L, 6L, 9L), 1L)
    if (L == 0L) return(NULL)
    data.frame(i = i, j = j, l = seq_len(L), d = rnorm(1), eta = rnorm(1))
  }))))
dat$y <- 2 * dat$d + dat$eta + rnorm(nrow(dat))

wl <- NULL
fl <- withCallingHandlers(
  with(dat, mwperm_layout(y = y, d = d, row = i, col = j, rep = l,
                          conf_int = FALSE, n_reps = 2L, seed = 1)),
  warning = function(w) {
    wl <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  })
ok(fl$pvalue == 1, grepl("no power", wl, fixed = TRUE),
   grepl("mwperm_irregular()", wl, fixed = TRUE))

fi2 <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                  L0 = 4L, min_block = 2L, conf_int = FALSE,
                                  n_reps = 2L, seed = 1))
ok(fi2$pvalue < 1, fi2$pvalue == 1 / (fi2$K + 1L),   # strongest attainable
   fi2$type == "irregular (Section 6.4)")
## no within-cell-variation warning here: constant d is the supported case
ok(is.null(withCallingHandlers({
  with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l, L0 = 4L,
                             min_block = 2L, conf_int = FALSE, n_reps = 1L,
                             seed = 1)); NULL },
  warning = function(w) stop("unexpected warning: ", conditionMessage(w)))))

## ---- 4. every retained cell holds exactly L0 observations ------------------
## Step (iii) of Section 6.4: the array handed to Procedure 2 is balanced.
ok(fi2$n_obs == fi2$cells_used * fi2$L0)
ok(any(grepl("Section 6.4", fi2$note, fixed = TRUE)))

## ---- 5. seeded reproducibility and dispatch equivalence --------------------
a <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                L0 = 4L, min_block = 2L, n_reps = 3L, seed = 5))
b <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                L0 = 4L, min_block = 2L, n_reps = 3L, seed = 5))
ok(identical(a$pvalue, b$pvalue), identical(a$conf_int, b$conf_int),
   identical(a$estimate, b$estimate))
u <- suppressWarnings(mwperm(y = "y", d = "d", index = c("i", "j"), rep = "l",
                             data = dat, design = "irregular", L0 = 4L,
                             min_block = 2L, n_reps = 3L, seed = 5,
                             verbose = FALSE))
ok(identical(u$pvalue, a$pvalue), identical(u$conf_int, a$conf_int),
   identical(u$auto$design, "irregular"))

## ---- 6. argument validation -------------------------------------------------
err(with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l)),
    "`L0` is required")
err(with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                               L0 = 1L)),
    "`L0` must be a single integer >= 2")
err(with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                               L0 = 500L)),
    "No cell has at least L0")
err(mwperm(y = "y", d = "d", index = c("i", "j"), rep = "l", data = dat,
           design = "irregular", verbose = FALSE),
    "requires `L0 =`")

cat("test-irregular.R: all assertions passed\n")
