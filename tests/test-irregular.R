## mwperm_irregular(): the Section 6.4 procedure.
##
## Section 6.4 as printed: (i) mask M_ij = 1{ell_ij >= L0}, Algorithm 2 on the
## mask, and in every retained cell drop ell_ij - L0 observations AT RANDOM;
## (ii) Procedure 2 on what remains, the within-cell position held fixed. The
## subsample is redrawn in every repetition and the reported p-value is the
## median (Remark 1; Appendix B repeats "the full permutation test" 100 times).
## That is `trim = "random"`, the default. `trim = "levels"` keeps instead a
## common set of L0 `rep` levels in every cell, which is what condition InvB
## needs when the repeats are periods with a common effect.
##
## The structural claim pinned here: on a complete balanced array with L0 equal
## to the number of within-cell slots nothing is dropped and Section 6.4
## degenerates to the panel construction of Section 6.2 -- the mask is all
## ones, the biclique search returns the whole array as one block, and
## Procedure 2's row/column groups act identically in every slot with the slot
## held fixed. This file pins the gather vectors themselves, not just a p-value.
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

## a group element as the engine applies it: the bare gather vector
bare <- function(op) { attributes(op) <- NULL; op }

## ---- 1. same permutations as the panel construction -----------------------
## The two front ends draw their groups from different sub-seed offsets (1, 2
## for the panel; 4q-1, 4q for the block builder), so the groups are matched by
## hand here and the resulting observation gather vectors compared directly.
## With L0 = T every cell already holds exactly L0 observations, so the random
## trim drops nothing and the position slot IS the period.
K <- m - 1L
rs <- 1L
Gr <- build_perm_set(m, K, seed = mwperm:::.sub_seed(rs, 3L))
Gc <- build_perm_set(n, K, seed = mwperm:::.sub_seed(rs, 4L))
A <- mwperm:::.build_obs_perms(cbind(ri, ci, ti), list(Gr, Gc, NULL),
                               design = "panel")
pdT <- with(tp, internal(".irregular_design")(row = importer, col = exporter,
                                              rep = year, L0 = TT,
                                              trim = "random",
                                              min_block = 2L,
                                              block_method = "greedy"))
opsT <- pdT$perm_builder(rs, K)
stopifnot(length(opsT) == K + 1L,
          identical(attr(opsT, "rows"), seq_len(nrow(tp))),   # nothing dropped
          identical(lapply(opsT, bare), A))

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
          fi$n_blocks == 1L, fi$cells_used == m * n, fi$L0 == TT,
          identical(fi$trim, "random"))

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

## ---- 4. the random trim: exactly L0 per retained cell, redrawn per rep ----
## Step (i) of Section 6.4 leaves every retained cell with exactly L0
## observations; each repetition draws its own subsample from its own seed, so
## the median over repetitions averages over the subsampling as Appendix B of
## the paper does.
stopifnot(fi2$n_obs == fi2$cells_used * fi2$L0,
          any(grepl("Section 6.4", fi2$note, fixed = TRUE)),
          any(grepl("at random", fi2$note, fixed = TRUE)),
          is.null(fi2$rep_levels))
pd4 <- with(dat, internal(".irregular_design")(row = i, col = j, rep = l,
                                               L0 = 4L, trim = "random",
                                               min_block = 2L,
                                               block_method = "greedy"))
cell_dat <- paste(dat$i, dat$j)
ell_dat <- ave(seq_len(nrow(dat)), cell_dat, FUN = length)
op1 <- pd4$perm_builder(1L, fi2$K)
op2 <- pd4$perm_builder(2L, fi2$K)
r1 <- attr(op1, "rows")                # rows of the engine data used by rep 1
r2 <- attr(op2, "rows")
stopifnot(length(r1) == fi2$n_obs, length(r2) == fi2$n_obs,
          !anyDuplicated(r1), identical(r1, sort(r1)),
          all(ell_dat[pd4$idx][r1] >= 4L),               # only cells that clear L0
          all(table(cell_dat[pd4$idx][r1]) == 4L),        # exactly L0 in each
          !identical(r1, r2),                             # redrawn per repetition
          identical(r1, attr(pd4$perm_builder(1L, fi2$K), "rows")))  # seeded
## every element is a bijection of the subsample that holds the within-cell
## position fixed and moves whole cells
pos1 <- ave(seq_along(r1), cell_dat[pd4$idx][r1], FUN = seq_along)
for (g in lapply(op1, bare)) {
  stopifnot(length(g) == length(r1), !anyDuplicated(g),
            identical(pos1[g], pos1))
}
for (g in lapply(op1[-1L], bare)) stopifnot(any(cell_dat[pd4$idx][r1][g] !=
                                                  cell_dat[pd4$idx][r1]))
## the seed = NULL path draws from the ambient stream and still balances
opN <- pd4$perm_builder(NULL, fi2$K)
stopifnot(all(table(cell_dat[pd4$idx][attr(opN, "rows")]) == 4L))

## ---- 5. seeded reproducibility, parallel = serial, dispatch equivalence ---
a <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                L0 = 4L, min_block = 2L, n_reps = 3L,
                                seed = 5))
b <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                L0 = 4L, min_block = 2L, n_reps = 3L,
                                seed = 5))
stopifnot(identical(a$pvalue, b$pvalue), identical(a$conf_int, b$conf_int),
          identical(a$estimate, b$estimate))
## the subsample is drawn inside the repetition loop, so the rep-parallel
## path must reproduce it from the rep seed
p2 <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                 L0 = 4L, min_block = 2L, n_reps = 3L,
                                 seed = 5, n_cores = 2L))
stopifnot(isTRUE(same_fit(a, p2)))
u <- suppressWarnings(mwperm(y = "y", d = "d", index = c("i", "j"), rep = "l",
                             data = dat, design = "irregular", L0 = 4L,
                             min_block = 2L, n_reps = 3L, seed = 5,
                             verbose = FALSE))
stopifnot(identical(u$pvalue, a$pvalue), identical(u$conf_int, a$conf_int),
          identical(u$auto$design, "irregular"))
## `trim` travels through mwperm() as well
uL <- suppressWarnings(mwperm(y = "y", d = "d", index = c("i", "j"), rep = "l",
                              data = dat, design = "irregular", L0 = 4L,
                              trim = "levels", min_block = 2L, n_reps = 3L,
                              seed = 5, verbose = FALSE))
aL <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                 L0 = 4L, trim = "levels", min_block = 2L,
                                 n_reps = 3L, seed = 5))
stopifnot(identical(uL$pvalue, aL$pvalue), identical(uL$conf_int, aL$conf_int),
          identical(aL$trim, "levels"))
## and `trim` is refused for the designs it does not apply to
expect_warn(mwperm(y = "log_trade", d = "log_dist",
                   index = c("importer", "exporter"), data = trade_dyadic,
                   trim = "levels", n_reps = 1L, seed = 1, conf_int = FALSE,
                   verbose = FALSE),
            "`trim` applies to the irregular design only")

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
expect_err(with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                      L0 = 4L, trim = "sometimes")),
           "'arg' should be one of")
expect_err(mwperm(y = "y", d = "d", index = c("i", "j"), rep = "l", data = dat,
                  design = "irregular", verbose = FALSE),
           "requires `L0 =`")
## the projection needs more retained observations per repetition than twice
## the nuisance columns; that count is cells_used x L0 (here 64 x 2 = 128),
## not the 256 observations of the retained cells (which would clear 2p = 130)
A0 <- expand.grid(t = 1:4, j = seq_len(8), i = seq_len(8))[, c("i", "j", "t")]
A0$d <- rnorm(nrow(A0)); A0$y <- rnorm(nrow(A0))
expect_err(with(A0, mwperm_irregular(y = y, d = d,
                                     x = matrix(rnorm(nrow(A0) * 64L),
                                                ncol = 64L),
                                     row = i, col = j, rep = t, L0 = 2L,
                                     min_block = 2L)),
           "Need N > 2p")

## ---- 7. trim = "levels": the slot is the rep LEVEL, held fixed across cells
## Condition InvB is (eps_ijt) =d (eps_[pi(i)][sigma(j)]t) with the SAME t on
## both sides. When the repeats are periods with a common effect, what the
## permutation holds fixed must be the `rep` level itself and every retained
## cell must carry the same set of levels; a random per-cell trim cannot
## guarantee that. The claim is pinned on the gather vectors: every element
## must map an observation to one with the SAME `rep`.
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
                               L0 = 2L, trim = "levels", min_block = 2L,
                               conf_int = FALSE, n_reps = 1L, seed = 1))
## only the cells observing BOTH retained levels survive: rows 1-6, all columns
stopifnot(fB$cells_used == (mB / 2) * nB, fB$n_obs == fB$cells_used * 2L,
          fB$cells_total == mB * nB, identical(fB$rep_levels, c("1", "2")),
          any(grepl("levels", fB$note, fixed = TRUE)))

## the front end's own retained set and permutation builder
pdB <- with(B, internal(".irregular_design")(row = i, col = j, rep = t, L0 = 2L,
                                             trim = "levels", min_block = 2L,
                                             block_method = "greedy"))
stopifnot(identical(pdB$S_labels, c("1", "2")),
          setequal(pdB$idx, which(B$i <= mB / 2)),
          length(pdB$idx) == fB$n_obs)
levB <- B$t[pdB$idx]
cellB <- paste(B$i, B$j)[pdB$idx]
opsB <- pdB$perm_builder(1L, fB$K)
stopifnot(length(opsB) == fB$K + 1L, identical(bare(opsB[[1L]]), seq_along(levB)),
          is.null(attr(opsB, "rows")))          # the cut is fixed, not per rep
for (g in lapply(opsB, bare)) stopifnot(identical(levB[g], levB))  # level kept
## ... and the check is not vacuous: every non-identity element moves cells
for (g in lapply(opsB[-1L], bare)) stopifnot(any(cellB[g] != cellB))

## The same design under the paper's random trim keeps every cell that has
## at least L0 observations (all 144 here: the count mask, not the level mask),
## and the slot it holds fixed is the position among the survivors -- which is
## NOT the period once cells keep different periods. Pinned so that the two
## modes cannot be confused.
fBr <- with(B, mwperm_irregular(y = y, d = d, row = i, col = j, rep = t,
                                L0 = 2L, min_block = 2L, conf_int = FALSE,
                                n_reps = 1L, seed = 1))
stopifnot(fBr$cells_used == mB * nB, fBr$n_obs == mB * nB * 2L,
          is.null(fBr$rep_levels))

## Design A: a rectangular array observed in periods 1..4 with L0 = 2. Under
## trim = "levels" the retained levels are {1, 2} in EVERY cell and every
## gather vector preserves the level; under trim = "random" each cell keeps a
## random pair of periods and the gather vectors preserve the position.
A <- expand.grid(t = 1:4, j = seq_len(8), i = seq_len(8))[, c("i", "j", "t")]
pdA <- with(A, internal(".irregular_design")(row = i, col = j, rep = t, L0 = 2L,
                                             trim = "levels", min_block = 2L,
                                             block_method = "greedy"))
stopifnot(identical(pdA$S_labels, c("1", "2")),
          setequal(pdA$idx, which(A$t <= 2L)),
          length(pdA$idx) == 2L * 64L)
levA <- A$t[pdA$idx]
for (g in lapply(pdA$perm_builder(7L, 7L), bare)) stopifnot(identical(levA[g], levA))
pdAr <- with(A, internal(".irregular_design")(row = i, col = j, rep = t, L0 = 2L,
                                              trim = "random", min_block = 2L,
                                              block_method = "greedy"))
opAr <- pdAr$perm_builder(7L, 7L)
rAr <- attr(opAr, "rows")
stopifnot(length(pdAr$idx) == 4L * 64L,                 # all 256 rows retained
          length(rAr) == 2L * 64L,                      # 2 per cell per rep
          all(table(paste(A$i, A$j)[pdAr$idx][rAr]) == 2L),
          !all(A$t[pdAr$idx][rAr] <= 2L))               # not the first two
posAr <- ave(seq_along(rAr), paste(A$i, A$j)[pdAr$idx][rAr], FUN = seq_along)
for (g in lapply(opAr, bare)) stopifnot(identical(posAr[g], posAr))

## rep = NULL under trim = "levels": the level is the within-cell order of
## appearance, so the common set is {1, ..., L0}, the mask is exactly GTW's
## 1{ell_ij >= L0}, and the first L0 observations of every cell that clears
## it are kept.
pdN <- with(dat, internal(".irregular_design")(row = i, col = j, rep = NULL,
                                               L0 = 4L, trim = "levels",
                                               min_block = 2L,
                                               block_method = "greedy"))
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
                                 trim = "levels", min_block = 2L,
                                 conf_int = FALSE, seed = 1))
fL <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                 L0 = 4L, trim = "levels", min_block = 2L,
                                 conf_int = FALSE, seed = 1))
stopifnot(isTRUE(same_fit(fN, fL, skip = c("call", "note"))))
## under the random trim `rep` only orders the survivors, and rep = NULL means
## order of appearance, so the two calls are identical there too
fNr <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, L0 = 4L,
                                  min_block = 2L, conf_int = FALSE, seed = 1))
fLr <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                  L0 = 4L, min_block = 2L, conf_int = FALSE,
                                  seed = 1))
stopifnot(isTRUE(same_fit(fNr, fLr, skip = c("call", "note"))))

## two observations in one cell at the same rep level contradict a shared
## level index (trim = "levels" only; the random trim never keys on the label)
B2 <- rbind(B, B[1L, ])
expect_err(with(B2, mwperm_irregular(y = y, d = d, row = i, col = j, rep = t,
                                     L0 = 2L, trim = "levels",
                                     min_block = 2L)),
           "more than once")

## ---- 8. the engine honours a per-repetition row subset --------------------
## The random trim reaches the engine as `attr(op, "rows")` on the builder's
## output: that repetition's statistics are computed on those rows only. A
## builder that returns the dyadic group over a fixed subset must therefore
## reproduce mwperm_dyadic() on that subset exactly.
data(trade_dyadic)
td <- trade_dyadic
sub <- with(td, which(as.integer(factor(importer)) <= 30L &
                        as.integer(factor(exporter)) <= 30L))
tds <- td[sub, ]
direct <- with(tds, mwperm_dyadic(y = log_trade, d = log_dist,
                                  x = cbind(log_gdp_i, log_gdp_j),
                                  row = importer, col = exporter,
                                  K = 9L, n_reps = 2L, seed = 3,
                                  conf_int = FALSE))
rr <- as.integer(factor(tds$importer)); cc <- as.integer(factor(tds$exporter))
builder <- function(rep_seed) {
  Gr <- build_perm_set(30L, 9L, seed = mwperm:::.sub_seed(rep_seed, 1L))
  Gc <- build_perm_set(30L, 9L, seed = mwperm:::.sub_seed(rep_seed, 2L))
  op <- mwperm:::.build_obs_perms(cbind(rr, cc), list(Gr, Gc))
  attr(op, "rows") <- sub
  op
}
hooked <- mwperm:::.ipt_engine(td$log_trade, as.matrix(td$log_dist),
                               cbind(1, td$log_gdp_i, td$log_gdp_j),
                               builder, K = 9L, n_reps = 2L, seed = 3,
                               alpha = 0.05, conf_int = FALSE, beta_null = 0,
                               grid = NULL, type = "dyadic", d_names = "d",
                               n_clusters = c(row = 30L, col = 30L),
                               call = quote(f()))
stopifnot(identical(hooked$pvalues_rep, direct$pvalues_rep),
          identical(hooked$pvalue, direct$pvalue))

passed("test-irregular.R")
