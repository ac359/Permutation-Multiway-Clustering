## mwperm_irregular(): the Section 6.4 procedure, step (i) as printed.
##
## (i) mask M_ij = 1{ell_ij >= L0}, Algorithm 2 on the mask, and in every
## retained cell drop ell_ij - L0 observations AT RANDOM; (ii) Procedure 2 on
## what remains, the within-cell position held fixed. The subsample is redrawn
## in every repetition and the reported p-value is the median (Remark 1;
## Appendix B repeats "the full permutation test" 100 times; the package
## default is 500). The random trim is exact when the observations inside a
## cell are exchangeable replicates; when the within-cell index is a period
## with an effect shared across cells it can over-reject even with a
## cell-constant `d` (the first note says so on every fit), a `d` varying over
## the periods makes that far worse (its own note, section 7), and the
## level-aligned cut that 0.4.1 offered as
## `trim = "levels"` is now mwperm_panel_missing(L0 = ) -- its assertions
## moved to tests/test-panel-missing.R section 8.
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
          !any(c("trim", "rep_levels") %in% names(fi)),   # retired in 0.4.2
          identical(formals(mwperm_irregular)$n_reps, 500L),
          !"trim" %in% names(formals(mwperm_irregular)))

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
          !any(c("trim", "rep_levels") %in% names(fi2)))
pd4 <- with(dat, internal(".irregular_design")(row = i, col = j, rep = l,
                                               L0 = 4L, min_block = 2L,
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
## `trim` was removed in 0.4.2: mwperm() keeps the slot for one release and
## names the replacement; the direct call has no such argument at all
expect_err(mwperm(y = "y", d = "d", index = c("i", "j"), rep = "l", data = dat,
                  design = "irregular", L0 = 4L, trim = "levels",
                  min_block = 2L, n_reps = 3L, seed = 5, verbose = FALSE),
           "`trim` was removed in 0.4.2")
expect_err(mwperm(y = "y", d = "d", index = c("i", "j"), rep = "l", data = dat,
                  design = "irregular", L0 = 4L, trim = "random",
                  min_block = 2L, n_reps = 3L, seed = 5, verbose = FALSE),
           "mwperm_panel_missing(time = <rep>, L0 = <L0>)")
expect_err(with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                      L0 = 4L, trim = "levels",
                                      min_block = 2L)),
           "unused argument")

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

## ---- 7. when `d` varies within cells: the note, and what the trim keeps -----
## The random trim leaves each cell with L0 observations but says nothing
## about WHICH L0. That is exact when the observations inside a cell are
## exchangeable replicates. When the within-cell index is a period with an
## effect shared across cells the test can over-reject even with a
## cell-constant `d` (size 0.39 with staggered observation windows), so the
## first note states the condition and names mwperm_panel_missing() on EVERY
## fit; when `d` also varies over the periods -- staggered adoption -- it is
## far worse (size 1.000, against 0.040 for the level-aligned cut), and the fit
## says so in a second note. Only `d` varying within a retained cell is
## detectable, so that is its trigger: a note, not a warning, because the fit
## IS valid for exchangeable replicates.
##
## Design B: rows 1-6 are observed in periods {1, 2}, rows 7-12 in {2, 3}, and
## d = 1{t >= start_ij} varies within each cell.
set.seed(3)
mB <- 12L; nB <- 12L
B <- do.call(rbind, lapply(seq_len(mB), function(i)
  do.call(rbind, lapply(seq_len(nB), function(j)
    data.frame(i = i, j = j, t = if (i <= mB / 2) 1:2 else 2:3)))))
startB <- matrix(sample(1:4, mB * nB, TRUE), mB, nB)
B$d <- as.numeric(B$t >= startB[cbind(B$i, B$j)])
B$y <- rnorm(mB)[B$i] + rnorm(nB)[B$j] + c(0, 0, 20)[B$t] + rnorm(nrow(B))
stopifnot(any(tapply(B$d, paste(B$i, B$j), function(v) diff(range(v)) > 0)))
fBr <- with(B, mwperm_irregular(y = y, d = d, row = i, col = j, rep = t,
                                L0 = 2L, min_block = 2L, conf_int = FALSE,
                                n_reps = 1L, seed = 1))
## the count mask keeps every cell with >= L0 observations (all 144: not the
## level mask), the slot held fixed is the position among the survivors, and
## the note fires
stopifnot(fBr$cells_used == mB * nB, fBr$n_obs == mB * nB * 2L,
          !any(c("trim", "rep_levels") %in% names(fBr)),
          any(grepl("`d` varies within cells", fBr$note, fixed = TRUE)),
          any(grepl("mwperm_panel_missing(time = <rep>, L0 = <L0>)", fBr$note,
                    fixed = TRUE)),
          any(grepl("NOT valid", fBr$note, fixed = TRUE)))
## it is a note, not a warning
stopifnot(length(warns_of(
  with(B, mwperm_irregular(y = y, d = d, row = i, col = j, rep = t, L0 = 2L,
                           min_block = 2L, conf_int = FALSE, n_reps = 1L,
                           seed = 1)))) == 0L)
## absent with a cell-constant d (sections 3-5's fixture), where the first
## note still states the condition and the route
stopifnot(!any(grepl("varies within cells", fi2$note, fixed = TRUE)),
          grepl("exact only when the observations inside a cell are ",
                fi2$note[1], fixed = TRUE),
          grepl("even with a cell-constant `d`", fi2$note[1], fixed = TRUE),
          grepl("mwperm_panel_missing(time = <rep>, L0 = <L0>)", fi2$note[1],
                fixed = TRUE))
## the trigger is `d` inside the RETAINED cells: a d that varies only in
## cells the mask drops does not fire it
dat_v <- dat
thin <- ave(seq_len(nrow(dat)), dat$i, dat$j, FUN = length) < 4L
dat_v$d[thin] <- dat_v$d[thin] + seq_len(sum(thin))
stopifnot(any(tapply(dat_v$d, paste(dat_v$i, dat_v$j),
                     function(v) diff(range(v)) > 0)))
fv <- with(dat_v, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                   L0 = 4L, min_block = 2L, conf_int = FALSE,
                                   n_reps = 1L, seed = 1))
stopifnot(!any(grepl("varies within cells", fv$note, fixed = TRUE)),
          isTRUE(same_fit(fv, fi2, skip = c("call", "note", "n_reps",
                                            "pvalues_rep", "pvalue"))))

## Design A: a rectangular array observed in periods 1..4 with L0 = 2: each
## cell keeps a random pair of periods and the gather vectors preserve the
## position among the survivors, not the period.
A <- expand.grid(t = 1:4, j = seq_len(8), i = seq_len(8))[, c("i", "j", "t")]
pdAr <- with(A, internal(".irregular_design")(row = i, col = j, rep = t, L0 = 2L,
                                              min_block = 2L,
                                              block_method = "greedy"))
opAr <- pdAr$perm_builder(7L, 7L)
rAr <- attr(opAr, "rows")
stopifnot(length(pdAr$idx) == 4L * 64L,                 # all 256 rows retained
          length(rAr) == 2L * 64L,                      # 2 per cell per rep
          all(table(paste(A$i, A$j)[pdAr$idx][rAr]) == 2L),
          !all(A$t[pdAr$idx][rAr] <= 2L))               # not the first two
posAr <- ave(seq_along(rAr), paste(A$i, A$j)[pdAr$idx][rAr], FUN = seq_along)
for (g in lapply(opAr, bare)) stopifnot(identical(posAr[g], posAr))

## `rep` only orders the survivors, and rep = NULL means order of appearance,
## so the two calls are identical
fNr <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, L0 = 4L,
                                  min_block = 2L, conf_int = FALSE,
                                  n_reps = 2L, seed = 1))
fLr <- with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
                                  L0 = 4L, min_block = 2L, conf_int = FALSE,
                                  n_reps = 2L, seed = 1))
stopifnot(isTRUE(same_fit(fNr, fLr, skip = c("call", "note"))))

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
