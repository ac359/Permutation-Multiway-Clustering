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

## ---- 1. same permutations as the panel construction -----------------------
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
stopifnot(identical(A, B))

## ---- 2. same data retained, hence same OLS reference ----------------------
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
stopifnot(identical(unname(fp$estimate), unname(fi$estimate)),
          identical(unname(fp$se_naive), unname(fi$se_naive)),
          identical(fp$n_obs, fi$n_obs), identical(fp$K, fi$K),
          fi$n_blocks == 1L, fi$cells_used == m * n, fi$L0 == TT)

## ---- 3. it works where the within-cell test has none ----------------------
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
stopifnot(fl$pvalue == 1, grepl("no power", wl, fixed = TRUE),
          grepl("mwperm_irregular()", wl, fixed = TRUE))

fi2 <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                  L0 = 4L, min_block = 2L, conf_int = FALSE,
                                  n_reps = 2L, seed = 1))
stopifnot(fi2$pvalue < 1,
          fi2$pvalue == 1 / (fi2$K + 1L),        # strongest attainable
          fi2$type == "irregular (Section 6.4)")
## no within-cell-variation warning here: constant d is the supported case
stopifnot(length(warns_of(
  with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                             L0 = 4L, min_block = 2L, conf_int = FALSE,
                             n_reps = 1L, seed = 1)))) == 0L)

## ---- 4. every retained cell holds exactly L0 observations -----------------
## Step (iii) of Section 6.4: the array handed to Procedure 2 is balanced.
stopifnot(fi2$n_obs == fi2$cells_used * fi2$L0)
stopifnot(any(grepl("Section 6.4", fi2$note, fixed = TRUE)))

## ---- 5. seeded reproducibility and dispatch equivalence -------------------
a <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                L0 = 4L, min_block = 2L, n_reps = 3L,
                                seed = 5))
b <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                L0 = 4L, min_block = 2L, n_reps = 3L,
                                seed = 5))
stopifnot(identical(a$pvalue, b$pvalue), identical(a$conf_int, b$conf_int),
          identical(a$estimate, b$estimate))
u <- suppressWarnings(mwperm(y = "y", d = "d", index = c("i", "j"), rep = "l",
                             data = dat, design = "irregular", L0 = 4L,
                             min_block = 2L, n_reps = 3L, seed = 5,
                             verbose = FALSE))
stopifnot(identical(u$pvalue, a$pvalue), identical(u$conf_int, a$conf_int),
          identical(u$auto$design, "irregular"))

## ---- 6. argument validation -----------------------------------------------
expect_err(with(dat, mwperm_irregular(y = y, d = d, row = i, col = j,
                                      rep = l)),
           "`L0` is required")
expect_err(with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                      L0 = 1L)),
           "`L0` must be a single integer >= 2")
expect_err(with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                      L0 = 500L)),
           "No cell has at least L0")
expect_err(mwperm(y = "y", d = "d", index = c("i", "j"), rep = "l", data = dat,
                  design = "irregular", verbose = FALSE),
           "requires `L0 =`")

## ---- 7. the slot is the rep LEVEL, held fixed across cells -----------------
## Section 6.4 is Procedure 2 combined with condition InvB of Section 6.2, and
## InvB is (eps_ijt) =d (eps_[pi(i)][sigma(j)]t) with the SAME t on both sides.
## So what the permutation holds fixed must be the `rep` level itself, and every
## retained cell must carry the same set of levels; otherwise "slot 2" is one
## period in one cell and another period in the next, and a common period
## effect is scrambled by the permutation. The claim is pinned on the gather
## vectors: every element must map an observation to one with the SAME `rep`.
##
## Design B: rows 1-6 are observed in periods {1, 2}, rows 7-12 in {2, 3}, so
## the cells do not share a period set, and d = 1{t >= start_ij} varies within
## each cell. The common set of L0 = 2 levels observed by the most cells is
## {1, 2} (tied with {2, 3}; the lower levels win), and only rows 1-6 clear it.
set.seed(3)
mB <- 12L; nB <- 12L
B <- do.call(rbind, lapply(seq_len(mB), function(i)
  do.call(rbind, lapply(seq_len(nB), function(j)
    data.frame(i = i, j = j, t = if (i <= mB / 2) 1:2 else 2:3)))))
startB <- matrix(sample(1:4, mB * nB, TRUE), mB, nB)
B$d <- as.numeric(B$t >= startB[cbind(B$i, B$j)])
B$y <- rnorm(mB)[B$i] + rnorm(nB)[B$j] + c(0, 0, 20)[B$t] + rnorm(nrow(B))
stopifnot(any(tapply(B$d, paste(B$i, B$j), function(v) diff(range(v)) > 0)))
fB <- with(B, mwperm_irregular(y = y, d = d, row = i, col = j, rep = t,
                               L0 = 2L, min_block = 2L, conf_int = FALSE,
                               n_reps = 1L, seed = 1))
## only the cells observing BOTH retained levels survive: rows 1-6, all columns
stopifnot(fB$cells_used == (mB / 2) * nB, fB$n_obs == fB$cells_used * 2L,
          fB$cells_total == mB * nB,
          any(grepl("levels", fB$note, fixed = TRUE)))

## the front end's own retained set and permutation builder
pdB <- with(B, internal(".irregular_design")(row = i, col = j, rep = t, L0 = 2L,
                                             min_block = 2L,
                                             block_method = "greedy"))
stopifnot(identical(pdB$S_labels, c("1", "2")),
          setequal(pdB$idx, which(B$i <= mB / 2)),
          length(pdB$idx) == fB$n_obs)
levB <- B$t[pdB$idx]
cellB <- paste(B$i, B$j)[pdB$idx]
opsB <- pdB$perm_builder(1L, fB$K)
stopifnot(length(opsB) == fB$K + 1L, identical(opsB[[1L]], seq_along(levB)))
for (g in opsB) stopifnot(identical(levB[g], levB))   # rep label preserved
## ... and the check is not vacuous: every non-identity element moves cells
for (g in opsB[-1L]) stopifnot(any(cellB[g] != cellB))

## Design A: a rectangular array observed in periods 1..4 with L0 = 2, so the
## trim itself must not break alignment. The retained levels are {1, 2} in
## EVERY cell -- a deterministic cut, not a per-cell random draw -- and every
## gather vector preserves the level.
A <- expand.grid(t = 1:4, j = seq_len(8), i = seq_len(8))[, c("i", "j", "t")]
pdA <- with(A, internal(".irregular_design")(row = i, col = j, rep = t, L0 = 2L,
                                             min_block = 2L,
                                             block_method = "greedy"))
stopifnot(identical(pdA$S_labels, c("1", "2")),
          setequal(pdA$idx, which(A$t <= 2L)),
          length(pdA$idx) == 2L * 64L)
levA <- A$t[pdA$idx]
for (g in pdA$perm_builder(7L, 7L)) stopifnot(identical(levA[g], levA))

## rep = NULL degeneracy: the level is the within-cell order of appearance, so
## the common set is {1, ..., L0}, the mask is exactly GTW's 1{ell_ij >= L0},
## and the first L0 observations of every cell that clears it are kept.
pdN <- with(dat, internal(".irregular_design")(row = i, col = j, rep = NULL,
                                               L0 = 4L, min_block = 2L,
                                               block_method = "greedy"))
ell_dat <- ave(seq_len(nrow(dat)), dat$i, dat$j, FUN = length)
pos_dat <- ave(seq_len(nrow(dat)), dat$i, dat$j, FUN = seq_along)
stopifnot(identical(pdN$S, 1:4),
          identical(pdN$n_mask, sum(tapply(ell_dat, paste(dat$i, dat$j),
                                           `[`, 1L) >= 4L)),
          all(pos_dat[pdN$idx] <= 4L),           # the FIRST L0 of each cell
          all(ell_dat[pdN$idx] >= 4L),           # ... of cells that clear L0
          identical(pdN$idx, which(ell_dat >= 4L & pos_dat <= 4L &
                                     pdN$in_block[pdN$cell] > 0L)))
## and with rep = order of appearance the fit is identical to rep = NULL
## (the note differs by the words "(order of appearance within the cell)")
fN <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, L0 = 4L,
                                 min_block = 2L, conf_int = FALSE, seed = 1))
fL <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                 L0 = 4L, min_block = 2L, conf_int = FALSE,
                                 seed = 1))
stopifnot(isTRUE(same_fit(fN, fL, skip = c("call", "note"))))

## two observations in one cell at the same rep level contradict a shared index
B2 <- rbind(B, B[1L, ])
expect_err(with(B2, mwperm_irregular(y = y, d = d, row = i, col = j, rep = t,
                                     L0 = 2L, min_block = 2L)),
           "more than once")

passed("test-irregular.R")
