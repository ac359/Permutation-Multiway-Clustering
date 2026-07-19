## 02_power.R -- Power as the true coefficient moves away from the null.
##
## Confirms the second half of the validity/power story: holding the design
## fixed, the IPT's rejection rate rises monotonically toward 1 as the true
## beta grows, so the size control of 01_size.R is not bought with a dead test.
## Dyadic DGP (Section 7.1), n = 30, normal covariate, phi2_err = 0.15 (the
## low-error-correlation cell of 01_size, where the signal-to-noise ratio
## makes the approach to power 1 visible; phi2_err = 0.9 is far noisier and
## power climbs much more slowly, though still monotonically).
##
## Usage:  Rscript 02_power.R           # 1500 sims/cell (default)
## Output: out/02_power.txt (+ out/02_power_summary.rds); cache under ./cache.
## Runtime ~1 min at the default N.

source("mc_lib.R")
suppressMessages(library(mwperm))
N <- as.integer(Sys.getenv("MC_N", "1500"))
BATCH <- 500L
sink_both("out/02_power.txt")
cat("==== 02 Power vs. the true effect size ====\n")
cat(sprintf("mwperm %s | %d sims/cell | dyadic n=30, normal, phi2_err=0.15 | %s\n\n",
            as.character(packageVersion("mwperm")), N, format(Sys.time())))

betas <- c(0, 0.05, 0.1, 0.15, 0.2)   # spans the transition from size to ~1
rows <- list()
for (b in betas) {
  sim <- function(s, dgp_seed, fit_seed) {
    dat <- dgp_dyadic71(30L, "normal", phi2_err = 0.15, beta = b)
    fit <- mwperm_dyadic(dat$y, dat$d, dat$x, dat$row, dat$col,
                         conf_int = FALSE, seed = fit_seed)
    c(p = fit$pvalue, est = unname(fit$estimate))
  }
  m <- mc_cell(sprintf("power_dyadic_b%03d_p15_v1", round(100 * b)), N, sim,
               params = list(design = "dyadic", n = 30, phi2 = 0.15,
                             beta = b), batch = BATCH)
  rej <- mean(m[, "p"] <= 0.05)
  cat(sprintf("  beta = %.2f   reject@.05 = %s%% (MC SE %s%%)   mean est = %+.3f\n",
              b, fmt_pct(rej),
              fmt_pct(sqrt(rej * (1 - rej) / N)), mean(m[, "est"])))
  rows[[length(rows) + 1L]] <- data.frame(beta = b, reject05 = rej,
                                          mean_est = mean(m[, "est"]))
}
tab <- do.call(rbind, rows)
saveRDS(tab, "out/02_power_summary.rds")
cat("\n---- summary ----\n")
mono <- all(diff(tab$reject05) > -3 * sqrt(0.25 / N))   # monotone within MC error
cat(sprintf("  power monotone increasing in beta : %s\n",
            if (mono) "yes" else "NO (INVESTIGATE)"))
cat(sprintf("  size at beta = 0                  : %s%%\n",
            fmt_pct(tab$reject05[tab$beta == 0])))
cat(sprintf("  power at the largest beta = %.2f   : %s%%\n",
            max(betas), fmt_pct(tab$reject05[tab$beta == max(betas)])))
cat("\nfull table:\n")
print(tab, digits = 3, row.names = FALSE)
sink()
