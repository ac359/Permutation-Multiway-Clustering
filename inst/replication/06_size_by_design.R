## 06_size_by_design.R -- Type-I error for the four designs 01_size.R omits:
## replicated layouts (Section 6.3), irregular layouts (Section 6.4),
## incomplete arrays (Procedure 2) and incomplete panels (InvB blockwise).
##
## Every cell runs at n_reps = 1 -- the configuration Theorem 1 covers -- and
## reports K next to the rejection rate, because a design whose floor 1/(K+1)
## exceeds alpha CANNOT reject and its measured size is a vacuous zero. Each
## cell asserts K + 1 >= 1/alpha before its size is quoted.
##
## The irregular design is run twice. The first arm is Section 6.4 as printed
## (the 0.4.1 default, trim = "random"): exchangeable replicates inside each
## cell, a covariate CONSTANT within each cell, a random per-cell trim to L0
## redrawn per repetition. The second arm has a covariate that VARIES within a
## cell (staggered adoption) and a common period effect in the errors: there
## the test is valid only if the permutation holds the `rep` LEVEL fixed in
## every cell, trim = "levels" -- under the random trim this arm rejects a true
## null essentially always (1.000 in 1000 simulations), which is the reason
## the argument exists.
##
## Usage:  Rscript 06_size_by_design.R          # 1000 sims/cell (default)
## Output: out/06_size_by_design.txt (+ _summary.rds); cache under ./cache.
## Runtime ~3-5 min at the default N.

source("mc_lib.R")
suppressMessages(library(mwperm))
N <- as.integer(Sys.getenv("MC_N", "1000"))
BATCH <- 250L
ALPHA <- 0.05
sink_both("out/06_size_by_design.txt")
cat("==== 06 Type-I error by design: layout, irregular, missing, panel ",
    "(bicliques) ====\n", sep = "")
cat(sprintf(paste0("mwperm %s | %d sims/cell | n_reps = 1 | rejection rule ",
                   "p <= %.2f | %s\n\n"),
            as.character(packageVersion("mwperm")), N, ALPHA,
            format(Sys.time())))

rows <- list()
report <- function(tag, m) {
  p <- m[, "p"]; K <- m[, "K"]
  Kmin <- min(K)
  attainable <- (Kmin + 1L) * ALPHA >= 1          # floor 1/(K+1) <= alpha
  st <- size_table(p, alphas = ALPHA)
  su <- if (length(unique(K)) == 1L) superuniformity(p, K) else NULL
  cat(sprintf(paste0("  %-34s K = %2d (floor %.4f, %s)  size@.05 = %s%% ",
                     "(CP upper %s%%)%s\n"),
              tag, Kmin, 1 / (Kmin + 1),
              if (attainable) "attainable" else "NOT ATTAINABLE",
              fmt_pct(st$size), fmt_pct(st$cp_upper1),
              if (is.null(su)) "" else paste0("  |  super-unif: ",
                if (attr(su, "any_violation")) "** VIOLATION **" else "clean")))
  if (!attainable)
    cat("    ^ the floor exceeds alpha, so this size certifies nothing\n")
  rows[[length(rows) + 1L]] <<- data.frame(
    cell = tag, K = Kmin, attainable = attainable, size05 = st$size,
    cp_lo = st$cp_lo, cp_upper1 = st$cp_upper1,
    su_violation = if (is.null(su)) NA else attr(su, "any_violation"))
}

## ---- layout (Section 6.3): within-cell permutation, L0 balancing -----------
## 10 x 10 cells holding 22-30 exchangeable replicates each; L0 = 22 balances
## them so K = 21. Errors: an arbitrary cell effect eta_ij plus i.i.d. noise,
## which is exchangeable within every cell. (A replicate effect zeta_l SHARED
## across cells would not be: independent per-cell permutations change its
## alignment across cells, so it lies outside what Section 6.3 permits -- see
## the README's layout paragraph.)
cat("---- replicated two-way layout, mwperm_layout(L0 = 22) ----\n")
sim_lay <- function(s, dgp_seed, fit_seed) {
  set.seed(dgp_seed)
  cells <- expand.grid(r = seq_len(10L), c = seq_len(10L))
  reps <- sample(22:30, nrow(cells), replace = TRUE)
  row <- rep(cells$r, reps); col <- rep(cells$c, reps)
  l <- unlist(lapply(reps, seq_len))
  cid <- rep(seq_len(nrow(cells)), reps)
  d <- rnorm(length(row))
  y <- rnorm(nrow(cells))[cid] + rnorm(length(row))
  fit <- mwperm_layout(y, d, row = row, col = col, rep = l, L0 = 22L,
                       n_reps = 1, seed = fit_seed, conf_int = FALSE)
  c(p = fit$pvalue, K = fit$K)
}
report("layout 10x10, 22-30 reps, L0 = 22",
       mc_cell("size06_layout_v2", N, sim_lay,
               params = list(design = "layout", L0 = 22L), batch = BATCH))

## ---- irregular (Section 6.4), covariate CONSTANT within cells --------------
## The software paper's Section 6.5 design: 30 x 30 cells, rows 26-30 thin
## (1-3 observations), the rest 4-9; L0 = 4 keeps the 25 x 30 block, K = 24.
cat("\n---- irregular layout, mwperm_irregular(L0 = 4), d constant in cells\n")
sim_irr_c <- function(s, dgp_seed, fit_seed) {
  set.seed(dgp_seed)
  cells <- expand.grid(i = seq_len(30L), j = seq_len(30L))
  thin <- cells$i > 25L
  ell <- ifelse(thin, sample(1:3, nrow(cells), replace = TRUE),
                sample(4:9, nrow(cells), replace = TRUE))
  row <- rep(cells$i, ell); col <- rep(cells$j, ell)
  l <- unlist(lapply(ell, seq_len)); cid <- rep(seq_len(nrow(cells)), ell)
  d <- rbinom(nrow(cells), 1, 0.4)[cid]
  y <- rnorm(nrow(cells))[cid] + rnorm(length(row))
  fit <- mwperm_irregular(y, d, row = row, col = col, rep = l, L0 = 4L,
                          min_block = 3L, n_reps = 1, seed = fit_seed,
                          conf_int = FALSE)
  c(p = fit$pvalue, K = fit$K)
}
report("irregular 30x30, cell-constant d",
       mc_cell("size06_irregular_const_v2", N, sim_irr_c,
               params = list(design = "irregular", L0 = 4L), batch = BATCH))

## ---- irregular (Section 6.4), covariate VARYING within cells ---------------
## 22 x 22 cells each observed in periods 1..4, L0 = 2 (so two periods are cut
## from every cell), d_ijt = 1{t >= start_ij} (staggered adoption), and a
## common period effect zeta = (0, 0, 0, 20) in the errors. Valid only when
## the retained periods are the SAME two in every cell: trim = "levels".
cat("\n---- irregular layout, mwperm_irregular(L0 = 2), d varies in cells\n")
sim_irr_w <- function(s, dgp_seed, fit_seed) {
  set.seed(dgp_seed)
  m <- 22L; TT <- 4L
  g <- expand.grid(t = seq_len(TT), j = seq_len(m), i = seq_len(m))
  start <- matrix(sample(seq_len(TT + 1L), m * m, TRUE), m, m)
  d <- as.numeric(g$t >= start[cbind(g$i, g$j)])
  y <- rnorm(m)[g$i] + rnorm(m)[g$j] + c(0, 0, 0, 20)[g$t] + rnorm(nrow(g))
  fit <- mwperm_irregular(y, d, row = g$i, col = g$j, rep = g$t, L0 = 2L,
                          trim = "levels", min_block = 20L, n_reps = 1,
                          seed = fit_seed, conf_int = FALSE)
  c(p = fit$pvalue, K = fit$K)
}
report("irregular 22x22x4, within-cell d, L0 = 2",
       mc_cell("size06_irregular_within_v2", N, sim_irr_w,
               params = list(design = "irregular", L0 = 2L), batch = BATCH))

## ---- missing (Procedure 2): 40 x 40 array with the diagonal deleted --------
## The no-self-trade mask of a gravity dataset (the software paper's Section
## 6.6 design): fixed, hence independent of the errors (Assumption 4), and the
## greedy search cuts it into two fully observed 20 x 20 blocks, so K = 19 and
## rejection at 0.05 is attainable. Two-way random-effect errors and a
## random-feature covariate, both from Eq. (13).
cat("\n---- incomplete array (no diagonal), mwperm_missing(min_block = 20)\n")
sim_mis <- function(s, dgp_seed, fit_seed) {
  set.seed(dgp_seed)
  n <- 40L
  ii <- rep(seq_len(n), times = n); jj <- rep(seq_len(n), each = n)
  d <- as.vector(rf_array(n, n, 0.4, 0.4))
  y <- 0.5 + as.vector(rf_array(n, n, 0.3, 0.3))
  keep <- ii != jj
  fit <- mwperm_missing(y[keep], d[keep], row = ii[keep], col = jj[keep],
                        min_block = 20L, n_reps = 1, seed = fit_seed,
                        conf_int = FALSE)
  c(p = fit$pvalue, K = fit$K)
}
report("missing 40x40, diagonal deleted",
       mc_cell("size06_missing_v2", N, sim_mis,
               params = list(design = "missing", n = 40L, mask = "no diag"),
               batch = BATCH))

## ---- incomplete panel (InvB blockwise): 26 x 26 x 4, six pairs thinned -----
## The panel DGP of mc_lib.R (arbitrary common trend, exchangeable across
## pairs within a period), with six random pairs observed in period 1 only.
## Those pairs fail the "observed in every period" mask, and the biclique
## search must drop at most six rows to recover a fully observed block, so a
## side of at least 20 -- K >= 19 -- is guaranteed.
cat("\n---- incomplete panel, mwperm_panel_missing(min_block = 20) ----\n")
sim_pm <- function(s, dgp_seed, fit_seed) {
  set.seed(dgp_seed)
  dat <- dgp_panel(26L, 4L, beta = 0, trend = TRUE)
  pair <- (dat$row - 1L) * 26L + dat$col
  thin <- sample(unique(pair), 6L)
  keep <- !(pair %in% thin & dat$time > 1L)
  fit <- mwperm_panel_missing(dat$y[keep], dat$d[keep], dat$x[keep, ],
                              row = dat$row[keep], col = dat$col[keep],
                              time = dat$time[keep], min_block = 20L,
                              n_reps = 1, seed = fit_seed, conf_int = FALSE)
  c(p = fit$pvalue, K = fit$K)
}
report("panel 26x26x4, 6 pairs thinned",
       mc_cell("size06_panel_missing_v1", N, sim_pm,
               params = list(design = "panel_missing", n = 26L, Tt = 4L),
               batch = BATCH))

tab <- do.call(rbind, rows)
saveRDS(tab, "out/06_size_by_design_summary.rds")
cat("\n---- summary ----\n")
## At n_reps = 1 the sizes sit near nominal rather than near 1%, so the
## one-sided CP upper bound of a 1000-sim cell is about size + 1.3 points and
## would flag a correct 4.8% as "> 5.5%". Evidence of OVER-rejection is what
## matters: a two-sided 95% CP interval lying entirely above alpha, or a
## super-uniformity violation anywhere on the p-value grid.
cat(sprintf("  cells where rejection at .05 is attainable       : %d / %d\n",
            sum(tab$attainable), nrow(tab)))
cat(sprintf("  attainable cells with evidence of over-rejection : %d / %d\n",
            sum(tab$attainable & tab$cp_lo > ALPHA), sum(tab$attainable)))
cat(sprintf("  cells with a super-uniformity violation          : %d / %d\n",
            sum(tab$su_violation, na.rm = TRUE),
            sum(!is.na(tab$su_violation))))
cat("\nfull table:\n")
print(tab, digits = 3, row.names = FALSE)
sink()
