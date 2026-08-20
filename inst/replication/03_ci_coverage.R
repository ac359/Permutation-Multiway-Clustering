## 03_ci_coverage.R -- Confidence-interval coverage.
##
## The confidence set is the test inverted, so its coverage is inherited from
## the test's validity: a nominal 95% interval should cover the true beta at
## least 95% of the time. Measured on the dyadic and panel designs at the true
## null beta = 0, with the package default n_reps = 10 (median aggregation).
##
## Usage:  Rscript 03_ci_coverage.R      # 1000 sims/cell (default)
## Output: out/03_ci_coverage.txt (+ _summary.rds); cache under ./cache.
## Runtime ~2-4 min at the default N (interval inversion is heavier per fit).

source("mc_lib.R")
suppressMessages(library(mwperm))
N <- as.integer(Sys.getenv("MC_N", "600"))
BATCH <- 200L
sink_both("out/03_ci_coverage.txt")
cat("==== 03 Confidence-interval coverage (nominal 95%) ====\n")
cat(sprintf("mwperm %s | %d sims/cell | true beta = 0 | default n_reps=10 | %s\n\n",
            as.character(packageVersion("mwperm")), N, format(Sys.time())))

rows <- list()

## ---- dyadic (n = 40 so a 95% interval is attainable: 1/(K+1) <= .05) --------
sim_dy <- function(s, dgp_seed, fit_seed) {
  dat <- dgp_dyadic71(40L, "normal", phi2_err = 0.9, beta = 0)
  fit <- mwperm_dyadic(dat$y, dat$d, dat$x, dat$row, dat$col, seed = fit_seed)
  ci <- fit$conf_int
  c(cover = as.numeric(ci[1] <= 0 && 0 <= ci[2]), width = diff(ci))
}
md <- mc_cell("cover_dyadic_n40_v1", N, sim_dy,
              params = list(design = "dyadic", n = 40), batch = BATCH)
cov_d <- mean(md[, "cover"]); ci_d <- cp_interval(sum(md[, "cover"]), N)
cat(sprintf("  dyadic n=40   coverage = %s%% (CP [%s, %s])   mean width = %.2f\n",
            fmt_pct(cov_d), fmt_pct(ci_d["lo"]), fmt_pct(ci_d["hi"]),
            mean(md[, "width"])))
rows[[1L]] <- data.frame(design = "dyadic n=40", coverage = cov_d,
                         cp_lo = ci_d["lo"], mean_width = mean(md[, "width"]),
                         row.names = NULL)

## ---- panel (n = 22 x T = 6, trending) ---------------------------------------
sim_pan <- function(s, dgp_seed, fit_seed) {
  dat <- dgp_panel(22L, 6L, beta = 0, trend = TRUE)
  fit <- mwperm_panel(dat$y, dat$d, dat$x, dat$row, dat$col, dat$time,
                      seed = fit_seed)
  ci <- fit$conf_int
  c(cover = as.numeric(ci[1] <= 0 && 0 <= ci[2]), width = diff(ci))
}
mp <- mc_cell("cover_panel_n22T6_v1", N, sim_pan,
              params = list(design = "panel", n = 22, Tt = 6), batch = BATCH)
cov_p <- mean(mp[, "cover"]); ci_p <- cp_interval(sum(mp[, "cover"]), N)
cat(sprintf("  panel 22xT6   coverage = %s%% (CP [%s, %s])   mean width = %.2f\n",
            fmt_pct(cov_p), fmt_pct(ci_p["lo"]), fmt_pct(ci_p["hi"]),
            mean(mp[, "width"])))
rows[[2L]] <- data.frame(design = "panel 22xT6", coverage = cov_p,
                         cp_lo = ci_p["lo"], mean_width = mean(mp[, "width"]),
                         row.names = NULL)

tab <- do.call(rbind, rows)
saveRDS(tab, "out/03_ci_coverage_summary.rds")
cat("\n---- summary ----\n")
cat(sprintf(paste0("  cells with coverage significantly above nominal",
                   " (CP lower > 95%%) : %d / %d\n"),
            sum(tab$cp_lo > 0.95), nrow(tab)))
cat("  (coverage at or above nominal is the valid direction)\n\nfull table:\n")
print(tab, digits = 3, row.names = FALSE)
sink()
