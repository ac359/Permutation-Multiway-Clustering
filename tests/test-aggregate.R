## Cross-repetition aggregation: one rule for the p-value and the interval.
##
## Theorem 1 gives finite-sample validity for a SINGLE random permutation group.
## The median of n_reps dependent randomised p-values (Remark 1) is a
## de-randomisation heuristic, not itself a level-alpha p-value; min(1, 2 x
## median) is. `aggregate` selects between them, and the SAME choice must drive
## the reported p-value, the confidence interval, and the joint region -- so no
## value the test accepts can fall outside the reported set. The choice also
## sets the REJECTION FLOOR: median2 can never return below 2/(K+1), which is
## what section 8 pins.
##
## Reaches internals via ::: -- run against a FRESHLY INSTALLED package.
library(mwperm)

ok <- function(...) stopifnot(...)

data(trade_dyadic)
td <- trade_dyadic

fit <- function(...) with(td, mwperm_dyadic(
  y = log_trade, d = log_dist, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, ...))

## ---- 1. median2 is exactly twice the median, capped at 1 -------------------
for (nr in c(1L, 5L, 10L)) {
  a <- fit(n_reps = nr, conf_int = FALSE, seed = 1)
  b <- fit(n_reps = nr, conf_int = FALSE, seed = 1, aggregate = "median2")
  ok(identical(b$pvalue, min(1, 2 * a$pvalue)))
  ok(identical(a$pvalues_rep, b$pvalues_rep))   # same permutations, same reps
}
## and it really does cap: a placebo p near 1 stays at 1
p <- with(td, mwperm_dyadic(y = log_trade, d = placebo,
                            x = cbind(log_gdp_i, log_gdp_j),
                            row = importer, col = exporter, conf_int = FALSE,
                            seed = 1, aggregate = "median2"))
ok(p$pvalue == 1)

## ---- 2. the default is unchanged -------------------------------------------
d1 <- fit(n_reps = 10L, seed = 1)
d2 <- fit(n_reps = 10L, seed = 1, aggregate = "median")
ok(identical(d1$pvalue, d2$pvalue), identical(d1$conf_int, d2$conf_int),
   identical(d1$conf_set, d2$conf_set))

## ---- 3. the median2 confidence set is never narrower ------------------------
m2 <- fit(n_reps = 10L, seed = 1, aggregate = "median2")
ok(m2$conf_int[1L] <= d1$conf_int[1L], m2$conf_int[2L] >= d1$conf_int[2L])

## ---- 4. the aggregation rule itself -----------------------------------------
P <- matrix(c(0.1, 0.2, 0.3,
              0.4, 0.5, 0.9), nrow = 2L, byrow = TRUE)
ok(identical(mwperm:::.agg_pvals(P, "median"), c(0.2, 0.5)),
   identical(mwperm:::.agg_pvals(P, "median2"), c(0.4, 1.0)),
   identical(mwperm:::.agg_pvals(P, "union"), c(0.3, 0.9)))
## a single rep is a no-op for median, still doubled for median2
P1 <- matrix(c(0.25, 0.7), ncol = 1L)
ok(identical(mwperm:::.agg_pvals(P1, "median"), c(0.25, 0.7)),
   identical(mwperm:::.agg_pvals(P1, "median2"), c(0.5, 1.0)))

## ---- 5. p-value and confidence set stay consistent under both rules --------
## The rule is applied in exactly one place for each, so a value strictly inside
## the set must not be rejected and a value outside it must be. (Strictly
## inside: the reported end points are the closure -- see test-exact-ci.R.)
for (agg in c("median", "median2")) {
  f <- fit(n_reps = 5L, alpha = 0.05, seed = 3, aggregate = agg)
  inb <- mean(f$conf_int)
  out <- f$conf_int[2L] + 0.5 * diff(f$conf_int)
  fi <- fit(n_reps = 5L, alpha = 0.05, seed = 3, aggregate = agg,
            beta_null = inb, conf_int = FALSE)
  fo <- fit(n_reps = 5L, alpha = 0.05, seed = 3, aggregate = agg,
            beta_null = out, conf_int = FALSE)
  ok(fi$pvalue > 0.05, fo$pvalue <= 0.05)
}

## ---- 6. every front end accepts it ------------------------------------------
data(trade_panel)
tp <- trade_panel
fp <- with(tp, mwperm_panel(y = log_trade, d = fta, row = importer,
                            col = exporter, time = year, conf_int = FALSE,
                            n_reps = 3L, seed = 1, aggregate = "median2"))
fpm <- with(tp, mwperm_panel(y = log_trade, d = fta, row = importer,
                             col = exporter, time = year, conf_int = FALSE,
                             n_reps = 3L, seed = 1))
ok(identical(fp$pvalue, min(1, 2 * fpm$pvalue)))
fu <- mwperm(y = "log_trade", d = "log_dist", x = c("log_gdp_i", "log_gdp_j"),
             index = c("importer", "exporter"), data = td, n_reps = 3L,
             seed = 1, aggregate = "median2", conf_int = FALSE,
             verbose = FALSE)
fum <- mwperm(y = "log_trade", d = "log_dist", x = c("log_gdp_i", "log_gdp_j"),
              index = c("importer", "exporter"), data = td, n_reps = 3L,
              seed = 1, conf_int = FALSE, verbose = FALSE)
ok(identical(fu$pvalue, min(1, 2 * fum$pvalue)))

## ---- 7. the joint region uses the same rule ---------------------------------
## Compared on a FIXED grid: the default region grid auto-expands, so a wider
## accepted set can be described by fewer (coarser) points and point counts are
## not comparable across runs. On identical candidates the comparison is
## exact: every point median accepts, median2 accepts too.
gg <- list(seq(-1.6, -0.2, length.out = 15L), seq(-0.6, 0.8, length.out = 15L))
r1 <- with(td, mwperm_dyadic(y = log_trade, d = cbind(log_dist, border),
                             x = cbind(log_gdp_i, log_gdp_j), row = importer,
                             col = exporter, n_reps = 3L, seed = 1, grid = gg))
r2 <- with(td, mwperm_dyadic(y = log_trade, d = cbind(log_dist, border),
                             x = cbind(log_gdp_i, log_gdp_j), row = importer,
                             col = exporter, n_reps = 3L, seed = 1, grid = gg,
                             aggregate = "median2"))
key <- function(m) apply(m, 1L, paste, collapse = "|")
ok(all(key(r1$conf_region) %in% key(r2$conf_region)),
   nrow(r2$conf_region) >= nrow(r1$conf_region))

## ---- 8. the rejection floor follows the aggregation rule --------------------
## median2 reports min(1, 2 * median), so its smallest attainable p-value is
## 2/(K+1), not the grid step 1/(K+1). At K = 19 and alpha = 0.05 the step is
## exactly alpha -- "median" can reject -- while "median2" can never return
## below 0.10, so gating the note and the confidence set on the step left the
## user with an unbounded interval and no explanation.
f_md <- fit(K = 19L, alpha = 0.05, n_reps = 3L, seed = 1)
f_m2 <- fit(K = 19L, alpha = 0.05, n_reps = 3L, seed = 1, aggregate = "median2")
ok(identical(f_md$resolution, f_m2$resolution),   # the GRID STEP is the same
   f_md$resolution == 1 / 20,
   f_md$p_floor == 1 / 20, f_m2$p_floor == 2 / 20)
## median: floor 0.05 <= alpha, so a set is computed and no resolution note
ok(!is.null(f_md$conf_int), all(is.finite(f_md$conf_int)),
   !any(grepl("smallest attainable", f_md$note)))
## median2: floor 0.10 > alpha, so no set, and the note names the right floor
## and the right requirement (K + 1 >= 2/alpha = 40)
ok(is.null(f_m2$conf_int),
   any(grepl("2/(K+1)", f_m2$note, fixed = TRUE)),
   any(grepl("K + 1 >= 40", f_m2$note, fixed = TRUE)),
   any(grepl("at least 40 levels", f_m2$note, fixed = TRUE)))
## confint() explains the floor rather than advising a pointless refit
msg <- tryCatch(confint(f_m2), error = conditionMessage)
ok(is.character(msg), grepl("2/(K+1)", msg, fixed = TRUE),
   grepl("at least 40 levels", msg, fixed = TRUE),
   !grepl("refit with conf_int = TRUE", msg, fixed = TRUE))
## with room under the floor, median2 gets its set back
f_m2b <- fit(K = 19L, alpha = 0.15, n_reps = 3L, seed = 1,
             aggregate = "median2")
ok(!is.null(f_m2b$conf_int), all(is.finite(f_m2b$conf_int)),
   !any(grepl("smallest attainable", f_m2b$note)))

cat("test-aggregate.R: all assertions passed\n")
