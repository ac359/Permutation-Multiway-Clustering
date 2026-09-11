## mwperm_panel_missing(): incomplete panels, condition InvB applied blockwise.
##
## The design claim is structural, and section 1 pins it on the gather vectors
## rather than on a p-value: given a COMPLETE array, this procedure must
## degenerate exactly to mwperm_panel(). The mask is then all ones, the
## biclique search returns the whole array as one block, nothing is discarded,
## and Procedure 2's row/column groups act identically in every period with the
## period held fixed -- which is precisely what mwperm_panel() does with time.
## If those two constructions ever disagree on a complete array, one of them is
## wrong.
##
## What the procedure adds beyond mwperm_panel() is the mask: a pair observed
## in some periods but not all is dropped whole (section 2), because the "same
## (pi, sigma) in every period" map is only defined on pairs present in every
## period.
##
## Reaches internals via ::: -- run against a FRESHLY INSTALLED package.
## Cross-cutting contracts live elsewhere; see tests/README.md.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

data(trade_panel)
tp <- trade_panel
ri <- as.integer(factor(tp$importer))
ci <- as.integer(factor(tp$exporter))
ti <- as.integer(factor(tp$year))
m <- max(ri); n <- max(ci); TT <- max(ti)

## ---- 1. on a complete array it IS the panel construction ------------------
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
  ri = ri, ci = ci, blk = rep(1L, length(ri)), lrow = ri, lcol = ci,
  slot = ti)
stopifnot(identical(A, B))

## ... and the fit keeps every observation, so the OLS reference agrees exactly
fp <- with(tp, mwperm_panel(log_trade, fta, x = cbind(log_gdp_i, log_gdp_j),
                            row = importer, col = exporter, time = year,
                            conf_int = FALSE, n_reps = 2, seed = 1))
fm <- with(tp, mwperm_panel_missing(log_trade, fta,
                                    x = cbind(log_gdp_i, log_gdp_j),
                                    row = importer, col = exporter,
                                    time = year, min_block = 2,
                                    conf_int = FALSE, n_reps = 2, seed = 1))
stopifnot(identical(unname(fp$estimate), unname(fm$estimate)),
          identical(unname(fp$se_naive), unname(fm$se_naive)),
          identical(fp$n_obs, fm$n_obs), identical(fp$K, fm$K),
          fm$n_blocks == 1L, fm$cells_used == m * n,
          fm$cells_total == m * n,
          identical(fm$type, "panel (bicliques)"))

## ---- 2. the mask keeps only pairs present in EVERY period -----------------
## A pair observed in some periods but not all is dropped whole -- not
## partially, and not silently aligned against a different pair's periods.
pair <- (ri - 1L) * n + ci
thin_pairs <- sort(unique(pair))[1:30]
drop_obs <- pair %in% thin_pairs & ti > 1L        # keep only period 1 for them
tp2 <- tp[!drop_obs, ]
f2 <- with(tp2, mwperm_panel_missing(log_trade, fta,
                                     x = cbind(log_gdp_i, log_gdp_j),
                                     row = importer, col = exporter,
                                     time = year, min_block = 3,
                                     conf_int = FALSE, n_reps = 2, seed = 1))
stopifnot(f2$cells_total == m * n - 0L,           # every pair is still present
          f2$cells_used <= m * n - length(thin_pairs),
          f2$n_obs == f2$cells_used * TT)         # kept cells are complete
## every retained observation belongs to a pair with all TT periods
keep_pair <- table(pair[!drop_obs])
stopifnot(all(keep_pair[as.character(thin_pairs)] == 1L))

## ---- 3. structure and seeded reproducibility ------------------------------
stopifnot(isTRUE(same_fit(
  f2, with(tp2, mwperm_panel_missing(log_trade, fta,
                                     x = cbind(log_gdp_i, log_gdp_j),
                                     row = importer, col = exporter,
                                     time = year, min_block = 3,
                                     conf_int = FALSE, n_reps = 2, seed = 1)),
  skip = "call")))
stopifnot(identical(names(f2$n_clusters), c("row", "col", "time")),
          identical(unname(f2$n_clusters), c(m, n, TT)),
          f2$n_blocks >= 1L,
          any(grepl("Incomplete panel", f2$note)))

## ---- 4. InvB: an arbitrary common time trend changes nothing --------------
## With time_fe = TRUE the period dummies span every zeta_t, so adding one to
## the outcome shifts y inside col(X) and, by FWL, cannot move the fit. This is
## the claim the design exists for, and it must survive the biclique
## restriction unchanged.
zeta <- c(0, 3.1, -2.4, 7.0, 1.5, -4.2)[as.integer(factor(tp2$year))]
f2z <- with(tp2, mwperm_panel_missing(log_trade + zeta, fta,
                                      x = cbind(log_gdp_i, log_gdp_j),
                                      row = importer, col = exporter,
                                      time = year, min_block = 3,
                                      conf_int = FALSE, n_reps = 2, seed = 1))
stopifnot(identical(f2$pvalue, f2z$pvalue),
          identical(f2$pvalues_rep, f2z$pvalues_rep),
          max(abs(f2$estimate - f2z$estimate)) < 1e-10)
## without the dummies the trend is not projected out and the estimate moves
f2n <- with(tp2, mwperm_panel_missing(log_trade + zeta, fta,
                                      x = cbind(log_gdp_i, log_gdp_j),
                                      row = importer, col = exporter,
                                      time = year, min_block = 3,
                                      time_fe = FALSE, conf_int = FALSE,
                                      n_reps = 2, seed = 1))
stopifnot(abs(f2n$estimate - f2$estimate) > 1e-6)

## ---- 5. row-order invariance ----------------------------------------------
set.seed(9)
o <- sample(nrow(tp2))
f2o <- with(tp2[o, ], mwperm_panel_missing(log_trade, fta,
                                           x = cbind(log_gdp_i, log_gdp_j),
                                           row = importer, col = exporter,
                                           time = year, min_block = 3,
                                           conf_int = FALSE, n_reps = 2,
                                           seed = 1))
stopifnot(identical(f2$pvalue, f2o$pvalue),
          abs(f2$estimate - f2o$estimate) < 1e-10)

## ---- 6. parallel == serial ------------------------------------------------
stopifnot(isTRUE(same_fit(
  f2, with(tp2, mwperm_panel_missing(log_trade, fta,
                                     x = cbind(log_gdp_i, log_gdp_j),
                                     row = importer, col = exporter,
                                     time = year, min_block = 3,
                                     conf_int = FALSE, n_reps = 2, seed = 1,
                                     n_cores = 2L)),
  skip = "call")))

## ---- 7. validation ---------------------------------------------------------
## A repeated (row, col, period) is an error: two observations would claim the
## same slot and the gather vector would resolve the collision silently.
expect_err(with(rbind(tp2, tp2[1L, ]),
                mwperm_panel_missing(log_trade, fta, row = importer,
                                     col = exporter, time = year)),
           "at most once")
## Nothing clears the mask: every pair is missing at least one period.
stagger <- tp[!(ti == (pair %% TT) + 1L), ]
expect_err(with(stagger,
                mwperm_panel_missing(log_trade, fta, row = importer,
                                     col = exporter, time = year)),
           "observed in all")
## Pairs clear the mask, but no block is big enough for the floor asked for.
expect_err(with(tp2,
                mwperm_panel_missing(log_trade, fta, row = importer,
                                     col = exporter, time = year,
                                     min_block = 40)),
           "min_block")
## the shared scalar-argument layer still applies
expect_err(with(tp2, mwperm_panel_missing(log_trade, fta, row = importer,
                                          col = exporter, time = year,
                                          min_block = 3, alpha = 0)),
           "`alpha`")

passed("test-panel-missing.R")
