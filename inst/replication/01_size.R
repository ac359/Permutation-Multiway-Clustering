## 01_size.R -- Type-I error under the null (the finite-sample-validity claim).
##
## Reproduces the package's headline result: the invariant permutation test
## (IPT) controls size at or below the nominal level in finite samples across
## every design, while the classical OLS t-test that ignores the multi-way
## clustering over-rejects grossly. Also confirms super-uniformity of the whole
## p-value distribution on the dyadic cells (validity at every attainable u).
##
## Usage (from this directory):
##   Rscript 01_size.R                 # 2000 sims/cell (default)
##   MC_N=10000 Rscript 01_size.R      # tighter Monte-Carlo error
## Output: out/01_size.txt (+ out/01_size_summary.rds); cache under ./cache
## (redirect with MWPERM_REPL_CACHE). Runtime ~1-2 min at the default N.

source("mc_lib.R")
suppressMessages(library(mwperm))
N <- as.integer(Sys.getenv("MC_N", "2000"))
BATCH <- 500L
sink_both("out/01_size.txt")
cat("==== 01 Type-I error under the null ====\n")
cat(sprintf("mwperm %s | %d sims/cell | rejection rule p <= alpha | %s\n\n",
            as.character(packageVersion("mwperm")), N, format(Sys.time())))

rows <- list()

## ---- dyadic: IPT vs the naive OLS t-test ------------------------------------
## n = 40 so the default K = n - 1 = 39 makes rejection at alpha = 0.05
## attainable at both 1/40 and 2/40 (a 95% interval is also feasible), giving
## a size near nominal rather than the sawtooth's most conservative point.
cat("---- dyadic (n = 40), IPT vs classical OLS t-test ----\n")
dy_cells <- list(
  list(cov = "normal",    phi2 = 0.15, tag = "normal/phi.15"),
  list(cov = "normal",    phi2 = 0.90, tag = "normal/phi.90"),
  list(cov = "lognormal", phi2 = 0.90, tag = "lognorm/phi.90"))
for (cc in dy_cells) {
  sim <- function(s, dgp_seed, fit_seed) {
    dat <- dgp_dyadic71(40L, cc$cov, phi2_err = cc$phi2, beta = 0)
    fit <- mwperm_dyadic(dat$y, dat$d, dat$x, dat$row, dat$col,
                         conf_int = FALSE, seed = fit_seed)
    c(p_ipt = fit$pvalue, K = fit$K,
      p_ols = naive_ols_p(dat$y, cbind(1, dat$x), dat$d))
  }
  m <- mc_cell(sprintf("size_dyadic_%s_phi%02d_n40_v1", cc$cov,
                       round(100 * cc$phi2)),
               N, sim, params = list(design = "dyadic", n = 40, cov = cc$cov,
                                     phi2 = cc$phi2), batch = BATCH)
  si <- size_table(m[, "p_ipt"]); so <- size_table(m[, "p_ols"])
  su <- superuniformity(m[, "p_ipt"], m[, "K"])
  i5 <- si$alpha == 0.05; o5 <- so$alpha == 0.05
  cat(sprintf(paste0("  %-15s IPT size@.05 = %s%% (CP upper %s%%)  |  ",
                     "OLS = %s%%  |  super-unif: %s\n"),
              cc$tag, fmt_pct(si$size[i5]), fmt_pct(si$cp_upper1[i5]),
              fmt_pct(so$size[o5]),
              if (attr(su, "any_violation")) "** VIOLATION **" else "clean"))
  rows[[length(rows) + 1L]] <- data.frame(
    design = "dyadic", cell = cc$tag, method = "IPT",
    size05 = si$size[i5], cp_upper1 = si$cp_upper1[i5],
    su_violation = attr(su, "any_violation"))
  rows[[length(rows) + 1L]] <- data.frame(
    design = "dyadic", cell = cc$tag, method = "naive OLS",
    size05 = so$size[o5], cp_upper1 = so$cp_upper1[o5], su_violation = NA)
}

## ---- panel and three-way: IPT controls size on its own valid DGP ------------
cat("\n---- panel (InvB) and three-way (InvA), IPT ----\n")
sim_pan <- function(s, dgp_seed, fit_seed) {
  dat <- dgp_panel(20L, 6L, beta = 0, trend = TRUE)
  fit <- mwperm_panel(dat$y, dat$d, dat$x, dat$row, dat$col, dat$time,
                      conf_int = FALSE, seed = fit_seed)
  c(p = fit$pvalue, K = fit$K)
}
mp <- mc_cell("size_panel_n20T6_v1", N, sim_pan,
              params = list(design = "panel", n = 20, Tt = 6), batch = BATCH)
sp <- size_table(mp[, "p"])
cat(sprintf("  %-15s IPT size@.05 = %s%% (CP upper %s%%)\n", "panel n20xT6",
            fmt_pct(sp$size[sp$alpha == 0.05]),
            fmt_pct(sp$cp_upper1[sp$alpha == 0.05])))
rows[[length(rows) + 1L]] <- data.frame(
  design = "panel", cell = "n20xT6 trend", method = "IPT",
  size05 = sp$size[sp$alpha == 0.05],
  cp_upper1 = sp$cp_upper1[sp$alpha == 0.05], su_violation = NA)

## 20x20x20 so the three-way default K = min(dims) - 1 = 19 makes rejection
## at alpha = 0.05 attainable (1/20 = 0.05); a smaller array can never reject
## at 0.05 and would report a trivial 0%.
sim_tw <- function(s, dgp_seed, fit_seed) {
  dat <- dgp_threeway(20L, 20L, beta = 0)
  fit <- mwperm_threeway(dat$y, dat$d, dat$x, dat$id1, dat$id2, dat$id3,
                         conf_int = FALSE, seed = fit_seed)
  c(p = fit$pvalue, K = fit$K)
}
mt <- mc_cell("size_threeway_n20T20_v1", N, sim_tw,
              params = list(design = "threeway", n = 20, Tt = 20),
              batch = BATCH)
stw <- size_table(mt[, "p"]); sutw <- superuniformity(mt[, "p"], mt[, "K"])
cat(sprintf("  %-15s IPT size@.05 = %s%% (CP upper %s%%)  |  super-unif: %s\n",
            "3-way 20^3", fmt_pct(stw$size[stw$alpha == 0.05]),
            fmt_pct(stw$cp_upper1[stw$alpha == 0.05]),
            if (attr(sutw, "any_violation")) "** VIOLATION **" else "clean"))
rows[[length(rows) + 1L]] <- data.frame(
  design = "threeway", cell = "20x20x20", method = "IPT",
  size05 = stw$size[stw$alpha == 0.05],
  cp_upper1 = stw$cp_upper1[stw$alpha == 0.05],
  su_violation = attr(sutw, "any_violation"))

tab <- do.call(rbind, rows)
saveRDS(tab, "out/01_size_summary.rds")
cat("\n---- summary ----\n")
cat(sprintf("  IPT cells with size@.05 > 5.5%% (CP upper)   : %d / %d\n",
            sum(tab$method == "IPT" & tab$cp_upper1 > 0.055),
            sum(tab$method == "IPT")))
cat(sprintf("  IPT cells with a super-uniformity violation : %d / %d\n",
            sum(tab$method == "IPT" & tab$su_violation, na.rm = TRUE),
            sum(tab$method == "IPT" & !is.na(tab$su_violation))))
cat(sprintf("  naive-OLS size range@.05                    : %s%% - %s%%\n",
            fmt_pct(min(tab$size05[tab$method == "naive OLS"])),
            fmt_pct(max(tab$size05[tab$method == "naive OLS"]))))
cat("\nfull table:\n")
print(tab, digits = 3, row.names = FALSE)
sink()
