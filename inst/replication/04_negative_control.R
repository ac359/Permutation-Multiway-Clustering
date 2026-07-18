## 04_negative_control.R -- Why panel is the default for 3-index arrays.
##
## A three-index array can be tested two ways. mwperm_threeway permutes all
## three dimensions (condition InvA), which is INVALID when the third dimension
## is time and the errors are autocorrelated over it -- permuting time destroys
## a real trend and the test over-rejects. mwperm_panel holds time fixed and
## permutes only the cross-section (condition InvB), which stays valid under an
## arbitrary common time trend. This script feeds the SAME trending-panel null
## to both tests: the three-way test should over-reject badly, the panel test
## should stay at or below nominal. That gap is the reason mwperm() dispatches
## to panel by default and only runs three-way on an explicit request.
##
## Usage:  Rscript 04_negative_control.R   # 1000 sims/cell (default)
## Output: out/04_negative_control.txt (+ _summary.rds); cache under ./cache.
## Runtime ~1 min at the default N.

source("mc_lib.R")
suppressMessages(library(mwperm))
N <- as.integer(Sys.getenv("MC_N", "1000"))
BATCH <- 500L
ALPHA <- 0.2                              # a lenient level makes the gap stark
sink_both("out/04_negative_control.txt")
cat("==== 04 Negative control: three-way vs panel on a trending panel ====\n")
cat(sprintf("mwperm %s | %d sims/cell | n=25 x T=6, common time trend, beta=0 | %s\n",
            as.character(packageVersion("mwperm")), N, format(Sys.time())))
cat(sprintf("rejection measured at alpha = %.2g\n\n", ALPHA))

sim <- function(s, dgp_seed, fit_seed) {
  dat <- dgp_panel(25L, 6L, beta = 0, trend = TRUE)
  ft <- mwperm_threeway(dat$y, dat$d, dat$x, dat$row, dat$col, dat$time,
                        conf_int = FALSE, seed = fit_seed)
  fp <- mwperm_panel(dat$y, dat$d, dat$x, dat$row, dat$col, dat$time,
                     conf_int = FALSE, seed = fit_seed)
  c(p_threeway = ft$pvalue, p_panel = fp$pvalue)
}
m <- mc_cell("negctrl_trendpanel_n25T6_v1", N, sim,
             params = list(n = 25, Tt = 6, trend = TRUE), batch = BATCH)
rej_tw <- mean(m[, "p_threeway"] <= ALPHA)
rej_pan <- mean(m[, "p_panel"] <= ALPHA)
cat(sprintf("  three-way (permutes time, INVALID here) : reject = %s%% (CP lower %s%%)\n",
            fmt_pct(rej_tw),
            fmt_pct(cp_lower1(sum(m[, "p_threeway"] <= ALPHA), N))))
cat(sprintf("  panel     (holds time, valid)           : reject = %s%% (CP upper %s%%)\n\n",
            fmt_pct(rej_pan),
            fmt_pct(cp_upper1(sum(m[, "p_panel"] <= ALPHA), N))))
tab <- data.frame(test = c("threeway", "panel"),
                  reject = c(rej_tw, rej_pan), alpha = ALPHA)
saveRDS(tab, "out/04_negative_control_summary.rds")
cat("---- summary ----\n")
cat(sprintf("  three-way over-rejects (reject >> alpha)      : %s\n",
            if (cp_lower1(sum(m[, "p_threeway"] <= ALPHA), N) > ALPHA)
              "yes -- invalid on a trend, as expected" else "no"))
cat(sprintf("  panel controls size (CP upper <= alpha + MC)  : %s\n",
            if (cp_upper1(sum(m[, "p_panel"] <= ALPHA), N) <= ALPHA + 0.02)
              "yes -- valid under an arbitrary trend" else "check"))
cat("\nfull table:\n")
print(tab, digits = 3, row.names = FALSE)
sink()
