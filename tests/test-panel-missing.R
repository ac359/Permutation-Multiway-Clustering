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

## ---- 8. L0: keep the best-covered L0 periods, the SAME ones in every cell --
## (0.4.2; the level-aligned cut that mwperm_irregular(trim = "levels") did
## in 0.4.1 now lives here, and these assertions moved with it.) Condition
## InvB is (eps_ijt) =d (eps_[pi(i)][sigma(j)]t) with the SAME t on both
## sides: what the permutation holds fixed must be the period itself, and
## every retained cell must carry the same set of periods. With L0 = NULL that
## set is every period; with an integer L0 it is the L0 periods jointly
## observed by the most pairs, chosen from the observation pattern alone
## (never from y -- Assumption 4), and the fit reports it in `periods_used`.
bare <- function(op) { attributes(op) <- NULL; op }
stopifnot(identical(names(formals(mwperm_panel_missing))[1:7],
                    c("y", "d", "x", "row", "col", "time", "L0")),
          is.null(formals(mwperm_panel_missing)$L0))

## L0 = n_t is the default mask, number for number (only the note and the
## reported period set differ)
f2L <- with(tp2, mwperm_panel_missing(log_trade, fta,
                                      x = cbind(log_gdp_i, log_gdp_j),
                                      row = importer, col = exporter,
                                      time = year, L0 = TT, min_block = 3,
                                      conf_int = FALSE, n_reps = 2, seed = 1))
stopifnot(isTRUE(same_fit(f2, f2L, skip = c("call", "note", "periods_used"))),
          is.null(f2$periods_used),
          identical(f2L$periods_used, as.character(sort(unique(tp2$year)))),
          any(grepl("L0 = 6", f2L$note, fixed = TRUE)))

## Design B: rows 1-6 are observed in periods {1, 2}, rows 7-12 in {2, 3}, so
## the pairs do not share a period set and NO pair clears the default mask;
## d = 1{t >= start_ij} varies within each cell. The common set of L0 = 2
## periods observed by the most pairs is {1, 2} (tied with {2, 3}; the lower
## periods win), and only rows 1-6 clear it.
set.seed(3)
mB <- 12L; nB <- 12L
B <- do.call(rbind, lapply(seq_len(mB), function(i)
  do.call(rbind, lapply(seq_len(nB), function(j)
    data.frame(i = i, j = j, t = if (i <= mB / 2) 1:2 else 2:3)))))
startB <- matrix(sample(1:4, mB * nB, TRUE), mB, nB)
B$d <- as.numeric(B$t >= startB[cbind(B$i, B$j)])
B$y <- rnorm(mB)[B$i] + rnorm(nB)[B$j] + c(0, 0, 20)[B$t] + rnorm(nrow(B))
stopifnot(any(tapply(B$d, paste(B$i, B$j), function(v) diff(range(v)) > 0)))
expect_err(with(B, mwperm_panel_missing(y, d, row = i, col = j, time = t,
                                        min_block = 2L)),
           "pass `L0 =`")                       # the way out is named
fB <- with(B, mwperm_panel_missing(y, d, row = i, col = j, time = t, L0 = 2L,
                                   min_block = 2L, time_fe = FALSE,
                                   conf_int = FALSE, n_reps = 1L, seed = 1))
## only the pairs observing BOTH retained periods survive: rows 1-6, all cols
stopifnot(fB$cells_used == (mB / 2) * nB, fB$n_obs == fB$cells_used * 2L,
          fB$cells_total == mB * nB, identical(fB$periods_used, c("1", "2")),
          identical(unname(fB$n_clusters[["time"]]), 2L),
          any(grepl("periods", fB$note, fixed = TRUE)))

## the front end's OWN retained set and permutation builder: every gather
## vector maps an observation to one with the SAME period
pdB <- with(B, internal(".panel_missing_design")(row = i, col = j, time = t,
                                                 L0 = 2L, min_block = 2L,
                                                 block_method = "greedy"))
stopifnot(identical(pdB$S_labels, c("1", "2")),
          setequal(pdB$idx, which(B$i <= mB / 2)),
          length(pdB$idx) == fB$n_obs)
levB <- B$t[pdB$idx]
cellB <- paste(B$i, B$j)[pdB$idx]
opsB <- pdB$perm_builder(1L, fB$K)
stopifnot(length(opsB) == fB$K + 1L,
          identical(bare(opsB[[1L]]), seq_along(levB)),
          is.null(attr(opsB, "rows")))          # the cut is fixed, not per rep
for (g in lapply(opsB, bare)) stopifnot(identical(levB[g], levB))  # period kept
## ... and the check is not vacuous: every non-identity element moves cells
for (g in lapply(opsB[-1L], bare)) stopifnot(any(cellB[g] != cellB))
## every retained cell observes ALL of S
stopifnot(all(tapply(levB, cellB, function(v) identical(sort(v), 1:2))))

## Design A: a rectangular 8 x 8 array observed in periods 1..4 with L0 = 2.
## The retained periods are {1, 2} in EVERY cell, every gather vector
## preserves the period -- and the fit is exactly the L0 = NULL fit on the
## data subset to those two periods, which is what "keep the best-covered
## periods" means.
A <- expand.grid(t = 1:4, j = seq_len(8), i = seq_len(8))[, c("i", "j", "t")]
set.seed(8)
A$d <- rnorm(nrow(A))
A$y <- rnorm(8)[A$i] + rnorm(8)[A$j] + c(0, 1, 5, 2)[A$t] + 0.3 * A$d +
  rnorm(nrow(A))
pdA <- with(A, internal(".panel_missing_design")(row = i, col = j, time = t,
                                                 L0 = 2L, min_block = 2L,
                                                 block_method = "greedy"))
stopifnot(identical(pdA$S_labels, c("1", "2")),
          setequal(pdA$idx, which(A$t <= 2L)),
          length(pdA$idx) == 2L * 64L)
levA <- A$t[pdA$idx]
for (g in lapply(pdA$perm_builder(7L, 7L), bare))
  stopifnot(identical(levA[g], levA))
fA2 <- with(A, mwperm_panel_missing(y, d, row = i, col = j, time = t, L0 = 2L,
                                    min_block = 2L, conf_int = FALSE,
                                    n_reps = 2L, seed = 1))
fAs <- with(A[A$t <= 2L, ],
            mwperm_panel_missing(y, d, row = i, col = j, time = t,
                                 min_block = 2L, conf_int = FALSE,
                                 n_reps = 2L, seed = 1))
stopifnot(isTRUE(same_fit(fA2, fAs, skip = c("call", "note", "periods_used"))),
          identical(fA2$periods_used, c("1", "2")), is.null(fAs$periods_used))

## seeded reproducibility and parallel == serial survive the cut
stopifnot(isTRUE(same_fit(fB, with(B, mwperm_panel_missing(
  y, d, row = i, col = j, time = t, L0 = 2L, min_block = 2L,
  time_fe = FALSE, conf_int = FALSE, n_reps = 1L, seed = 1)), skip = "call")))
fB3 <- with(B, mwperm_panel_missing(y, d, row = i, col = j, time = t, L0 = 2L,
                                    min_block = 2L, time_fe = FALSE,
                                    conf_int = FALSE, n_reps = 3L, seed = 1))
stopifnot(isTRUE(same_fit(fB3, with(B, mwperm_panel_missing(
  y, d, row = i, col = j, time = t, L0 = 2L, min_block = 2L,
  time_fe = FALSE, conf_int = FALSE, n_reps = 3L, seed = 1, n_cores = 2L)),
  skip = "call")))

## validation: L0 is an integer in [2, number of periods]
expect_err(with(B, mwperm_panel_missing(y, d, row = i, col = j, time = t,
                                        L0 = 1L, min_block = 2L)),
           "`L0` must be a single integer >= 2")
expect_err(with(B, mwperm_panel_missing(y, d, row = i, col = j, time = t,
                                        L0 = c(2L, 3L), min_block = 2L)),
           "`L0` must be a single integer >= 2")
expect_err(with(B, mwperm_panel_missing(y, d, row = i, col = j, time = t,
                                        L0 = 4L, min_block = 2L)),
           "`L0`")                              # only 3 periods exist
## L0 periods exist but no pair observes L0 of them
one_each <- B[B$t == ifelse(B$i <= mB / 2, 1L, 3L), ]
expect_err(with(one_each, mwperm_panel_missing(y, d, row = i, col = j,
                                               time = t, L0 = 2L,
                                               min_block = 2L)),
           "L0 = 2")

passed("test-panel-missing.R")
