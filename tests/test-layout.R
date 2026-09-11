## mwperm_layout() -- the replicated two-way layout, within-cell
## exchangeability.
##
## Each (i, j) cell holds ell_ij replicates and permutations act only WITHIN
## cells, drawn independently per cell (Guo, Toulis & Wang 2026, Section 6.3).
## The group order is therefore set by the smallest cell, not by the number of
## clusters, and `L0` (Section 6.4) buys resolution back by dropping thin cells
## and downsampling the rest to a common size.
##
## Two things here are easy to get silently wrong and are pinned hard: the L0
## downsample is seeded (an unseeded one would make results irreproducible),
## and a covariate that is constant within every cell has no within-cell
## variation left after residualizing, so the test has no power and must say so.
##
## Cross-cutting contracts live elsewhere; see tests/README.md.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

## ---- fixture: 4 x 4 cells, 5 replicates each ------------------------------
set.seed(3)
gl <- expand.grid(l = 1:5, i = 1:4, j = 1:4)
Nl <- nrow(gl)
cell <- (gl$i - 1L) * 4L + gl$j
dl <- rnorm(Nl)
yl <- rnorm(16)[cell] + 0.5 * dl + rnorm(Nl)      # arbitrary cell effects

## ---- 1. structure and the within-cell group order -------------------------
## K defaults to min(cell size) - 1: the permutation is over replicates inside
## a cell, so the smallest cell caps it. n_clusters reports the layout, not a
## number of clusters per dimension.
fit <- mwperm_layout(yl, dl, row = gl$i, col = gl$j, seed = 1,
                     conf_int = FALSE)
stopifnot(inherits(fit, "mwperm"),
          identical(fit$type, "layout"),
          fit$K == 4L, fit$n_perm == 5L, fit$n_obs == Nl,
          identical(names(fit$n_clusters), c("ncell", "min_cell")),
          identical(unname(fit$n_clusters), c(16L, 5L)))
stopifnot(isTRUE(same_fit(fit, mwperm_layout(yl, dl, row = gl$i, col = gl$j,
                                             seed = 1, conf_int = FALSE))))

## arbitrary cell effects are what the within-cell argument is for: the design
## carries a strong signal and the p-value reaches its floor
stopifnot(identical(
  mwperm_layout(yl + 3 * dl, dl, row = gl$i, col = gl$j, seed = 1,
                conf_int = FALSE)$pvalue,
  fit$p_floor))

## ---- 2. unequal cell sizes, and L0 balancing ------------------------------
## Without L0 the cells keep their own sizes and K follows the smallest. With
## L0, thin cells are dropped and the rest are downsampled to exactly L0, so
## K = L0 - 1 and the note says the array was balanced.
gu <- expand.grid(l = 1:6, i = 1:5, j = 1:5)
gu <- gu[!(gu$l > 3 & gu$i == 1), ]                # cells in row 1 hold 3
Nu <- nrow(gu)
set.seed(202)
yu <- rnorm(Nu)
du <- rnorm(Nu)
f_raw <- mwperm_layout(yu, du, row = gu$i, col = gu$j, seed = 2, n_reps = 2,
                       conf_int = FALSE)
stopifnot(f_raw$K == 2L,                           # smallest cell holds 3
          f_raw$n_clusters[["min_cell"]] == 3L)
f_L0 <- mwperm_layout(yu, du, row = gu$i, col = gu$j, L0 = 3, seed = 2,
                      n_reps = 2, conf_int = FALSE)
stopifnot(any(grepl("Balanced to L0", f_L0$note)), f_L0$K <= 2L)

## the downsample must be seeded: two identical calls agree exactly
f_L0b <- mwperm_layout(yu, du, row = gu$i, col = gu$j, L0 = 3, seed = 2,
                       n_reps = 2, conf_int = FALSE)
stopifnot(isTRUE(same_fit(f_L0, f_L0b)))
## raising L0 above every cell size leaves nothing to test, and says so
expect_err(mwperm_layout(yu, du, row = gu$i, col = gu$j, L0 = 99), "L0 = 99")
expect_err(mwperm_layout(yu, du, row = gu$i, col = gu$j, L0 = 1), "`L0`")

## ---- 3. `rep` labels the replicates, and makes the fit order-free ---------
## With `rep` supplied the within-cell slot of an observation is its label, not
## its position in the input, so shuffling the rows cannot move it.
f0 <- mwperm_layout(yl, dl, row = gl$i, col = gl$j, rep = gl$l, seed = 3,
                    conf_int = FALSE)
set.seed(77)
o <- sample(Nl)
f1 <- mwperm_layout(yl[o], dl[o], row = gl$i[o], col = gl$j[o],
                    rep = gl$l[o], seed = 3, conf_int = FALSE)
stopifnot(identical(f0$pvalue, f1$pvalue),
          abs(f0$estimate - f1$estimate) < 1e-10)
expect_err(mwperm_layout(yl, dl, row = gl$i, col = gl$j, rep = gl$l[-1]),
           "`rep`")
## an NA `rep` used to be ranked silently last within its cell
repNA <- as.numeric(gl$l)
repNA[3L] <- NA
expect_err(mwperm_layout(yl, dl, row = gl$i, col = gl$j, rep = repNA), "`rep`")

## ---- 4. a cell-level covariate has no within-cell variation ---------------
## Residualizing d on the cell means leaves nothing, so the test has no power
## whatever the effect size. That is a property of the design, not a failure,
## and the front end must warn rather than report a confident-looking p-value.
d_cell <- rnorm(25)[(gu$i - 1L) * 5L + gu$j]      # constant inside each cell
expect_warn(mwperm_layout(yu, d_cell, row = gu$i, col = gu$j, seed = 1,
                          n_reps = 1, conf_int = FALSE),
            "constant within every cell")

passed("test-layout.R")
