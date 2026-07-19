## 05_permute.R -- One-sided permutation for incomplete arrays (permute=).
##
## mwperm_missing(permute = "rows"/"cols") permutes a single dimension, a
## subgroup of the full invariance group, so the test stays valid under the
## same exchangeability assumption while K + 1 is capped by the permuted block
## side alone. This lets designs whose fully observed blocks are short in one
## dimension -- down to a single fully observed column -- reach resolutions the
## two-sided test cannot. This script checks size control and power for the
## one-sided test on exactly such a design.
##
## Usage:  Rscript 05_permute.R          # 2000 sims (size), 500 (power) default
## Output: out/05_permute.txt (+ _summary.rds); cache under ./cache.
## Runtime ~1 min at the default N.

source("mc_lib.R")
suppressMessages(library(mwperm))
N <- as.integer(Sys.getenv("MC_N", "2000"))
BATCH <- 500L
sink_both("out/05_permute.txt")
cat("==== 05 One-sided permutation (permute = 'rows') ====\n")
cat(sprintf("mwperm %s | size %d sims, power 500 sims | %s\n\n",
            as.character(packageVersion("mwperm")), N, format(Sys.time())))

## ---- size: a single fully observed column (30 rows x 1 column) --------------
## Two-sided permutation has NO usable block here; rows-only reaches K = 29
## (resolution 1/30). Errors i.i.d. across the exchangeable rows.
n1 <- 30L
sim_size <- function(s, dgp_seed, fit_seed) {
  set.seed(dgp_seed)
  d <- rnorm(n1); y <- 0.5 + rnorm(n1)                 # beta = 0 (null)
  fit <- mwperm_missing(y, d, row = seq_len(n1), col = rep(1L, n1),
                        permute = "rows", n_reps = 1, seed = fit_seed,
                        conf_int = FALSE)
  c(p = fit$pvalue, K = fit$K)
}
ms <- mc_cell("permute_rows_size_col1_v1", N, sim_size,
              params = list(n = n1), batch = BATCH)
ss <- size_table(ms[, "p"]); su <- superuniformity(ms[, "p"], ms[, "K"])
cat("---- size (30x1 single-column, rows-only, null) ----\n")
cat(sprintf("  K = %d  |  size@.05 = %s%% (CP upper %s%%)  |  super-unif: %s\n\n",
            ms[1, "K"], fmt_pct(ss$size[ss$alpha == 0.05]),
            fmt_pct(ss$cp_upper1[ss$alpha == 0.05]),
            if (attr(su, "any_violation")) "** VIOLATION **" else "clean"))

## ---- power: same design, beta away from 0 -----------------------------------
cat("---- power (30x1 single-column, rows-only) ----\n")
betas <- c(0.25, 0.5, 1.0)
pw <- vapply(betas, function(b) {
  sim <- function(s, dgp_seed, fit_seed) {
    set.seed(dgp_seed)
    d <- rnorm(n1); y <- 0.5 + b * d + rnorm(n1)
    fit <- mwperm_missing(y, d, row = seq_len(n1), col = rep(1L, n1),
                          permute = "rows", n_reps = 1, seed = fit_seed,
                          conf_int = FALSE)
    c(p = fit$pvalue)
  }
  m <- mc_cell(sprintf("permute_rows_power_b%03d_v1", round(100 * b)),
               500L, sim, params = list(n = n1, beta = b), batch = 500L)
  mean(m[, "p"] <= 0.05)
}, numeric(1))
pw0 <- mean(ms[, "p"] <= 0.05)
for (i in seq_along(betas))
  cat(sprintf("  beta = %.2f   reject@.05 = %s%%\n", betas[i], fmt_pct(pw[i])))

tab <- data.frame(beta = c(0, betas), reject05 = c(pw0, pw))
saveRDS(list(size = ss, power = tab), "out/05_permute_summary.rds")
cat("\n---- summary ----\n")
cat(sprintf("  size controlled (CP upper <= 5.5%%)        : %s\n",
            if (ss$cp_upper1[ss$alpha == 0.05] <= 0.055) "yes" else "check"))
cat(sprintf("  power monotone, ~1 at beta = 1            : %s\n",
            if (all(diff(c(pw0, pw)) > 0) && pw[length(pw)] > 0.9)
              "yes" else "check"))
cat("\npower table:\n")
print(tab, digits = 3, row.names = FALSE)
sink()
