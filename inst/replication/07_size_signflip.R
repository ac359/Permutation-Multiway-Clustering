## 07_size_signflip.R -- The sign-flip test (IPT-Het) against the permutation
## test under heteroskedastic errors.
##
## Every permutation design assumes the errors exchangeable GIVEN the
## covariates, which an error variance that depends on the covariates
## violates. This script runs mwperm_dyadic() and mwperm_dyadic_het() on the
## same draws of the paper authors' heteroskedastic gravity design (n = 25
## clusters per side; log GDP node effects; a two-way random-feature log
## distance; error sd = 0.35 * exp(rho1 * standardised gravity mean + rho2 *
## standardised distance); centred lognormal innovations, i.e. ASYMMETRIC) at
## three strengths of heteroskedasticity under the null, and at beta = 0.15
## under homoskedasticity for the power comparison.
##
## What it shows: the permutation test's rejection rate rises with the
## heteroskedasticity while the sign-flip test's stays at or below nominal,
## and under homoskedastic errors the sign-flip test is the less powerful of
## the two. It is NOT a strict upgrade -- that trade-off is the point.
##
## Both tests run at n_reps = 1 (Theorem 1's configuration): the permutation
## test at its default K = 24, the sign-flip test at n_flip = 6 (group order
## 32, so rejection at 0.05 is attainable; the default 8 costs four times as
## much and, at this n, changes little).
##
## Usage:  Rscript 07_size_signflip.R          # 500 sims/cell (default)
## Output: out/07_size_signflip.txt (+ _summary.rds); cache under ./cache.
## Runtime ~7 min at the default N (two fits per sim).

source("mc_lib.R")
suppressMessages(library(mwperm))
N <- as.integer(Sys.getenv("MC_N", "500"))
BATCH <- 250L
ALPHA <- 0.05
sink_both("out/07_size_signflip.txt")
cat("==== 07 Sign-flip test (IPT-Het) vs permutation test under ",
    "heteroskedasticity ====\n", sep = "")
cat(sprintf(paste0("mwperm %s | %d sims/cell | n = 25 | perm K = 24, flip ",
                   "n_flip = 6 | n_reps = 1 | rejection rule p <= %.2f | %s\n\n"),
            as.character(packageVersion("mwperm")), N, ALPHA,
            format(Sys.time())))

## The authors' simulate_gravity_model() in long format (i varies fastest).
dgp_gravity_het <- function(n, b, rho1, rho2, phi1 = 0.4, phi2 = 0.4) {
  m <- n
  log_gdp <- rnorm(n, mean = log(500), sd = 1.2)
  log_gdp <- log(pmin(pmax(exp(log_gdp), 5), 25000))
  s1sq <- phi1 / ((1 - phi1) * (1 - phi2) - phi1 * phi2)
  s2sq <- phi2 / ((1 - phi1) * (1 - phi2) - phi1 * phi2)
  vg <- rnorm(m); vh <- rnorm(n)
  g <- expand.grid(i = seq_len(m), j = seq_len(n))
  dist_raw <- sqrt(s1sq) * vg[g$i] + sqrt(s2sq) * vh[g$j] + rnorm(m * n)
  D <- log(3000) + 0.5 * as.numeric(scale(dist_raw))
  X <- cbind(log_gdp[g$i], log_gdp[g$j])
  mu <- 4 + 0.8 * X[, 1L] + 0.8 * X[, 2L] + b * D
  sigma <- 0.35 * exp(rho1 * as.numeric(scale(mu)) + rho2 * as.numeric(scale(D)))
  u <- rlnorm(m * n, meanlog = 0, sdlog = 1)
  eta <- (u - exp(0.5)) / sqrt((exp(1) - 1) * exp(1))
  list(y = mu + sigma * eta, d = D, x = X, row = g$i, col = g$j)
}

cells <- list(
  list(tag = "null, strong het (rho 0.9, 0.6)", rho1 = 0.9, rho2 = 0.6, b = 0),
  list(tag = "null, mild het   (rho 0.5, 0.3)", rho1 = 0.5, rho2 = 0.3, b = 0),
  list(tag = "null, homoskedastic",             rho1 = 0.0, rho2 = 0.0, b = 0),
  list(tag = "power, homoskedastic, beta = .15", rho1 = 0.0, rho2 = 0.0,
       b = 0.15))

rows <- list()
for (ci in seq_along(cells)) {
  cl <- cells[[ci]]
  sim <- function(s, dgp_seed, fit_seed) {
    dat <- dgp_gravity_het(25L, b = cl$b, rho1 = cl$rho1, rho2 = cl$rho2)
    fp <- mwperm_dyadic(dat$y, dat$d, x = dat$x, row = dat$row, col = dat$col,
                        n_reps = 1, seed = fit_seed, conf_int = FALSE)
    ff <- mwperm_dyadic_het(dat$y, dat$d, x = dat$x, row = dat$row,
                            col = dat$col, n_flip = 6L, n_reps = 1,
                            seed = fit_seed, conf_int = FALSE)
    c(p_perm = fp$pvalue, K_perm = fp$K, p_flip = ff$pvalue, K_flip = ff$K)
  }
  m <- mc_cell(sprintf("size07_signflip_cell%d_v1", ci), N, sim,
               params = list(rho1 = cl$rho1, rho2 = cl$rho2, b = cl$b,
                             n = 25L, n_flip = 6L), batch = BATCH)
  st_p <- size_table(m[, "p_perm"], alphas = ALPHA)
  st_f <- size_table(m[, "p_flip"], alphas = ALPHA)
  stopifnot((min(m[, "K_perm"]) + 1L) * ALPHA >= 1,     # floors attainable
            (min(m[, "K_flip"]) + 1L) * ALPHA >= 1)
  cat(sprintf(paste0("  %-34s perm (K = %2d): %s%%  [CP %s, %s]   ",
                     "flip (order %3d): %s%%  [CP %s, %s]\n"),
              cl$tag, min(m[, "K_perm"]), fmt_pct(st_p$size),
              fmt_pct(st_p$cp_lo), fmt_pct(st_p$cp_upper1),
              min(m[, "K_flip"]) + 1L, fmt_pct(st_f$size),
              fmt_pct(st_f$cp_lo), fmt_pct(st_f$cp_upper1)))
  rows[[ci]] <- data.frame(cell = cl$tag, rho1 = cl$rho1, rho2 = cl$rho2,
                           beta = cl$b, rej_perm = st_p$size,
                           rej_flip = st_f$size,
                           cp_lo_perm = st_p$cp_lo, cp_lo_flip = st_f$cp_lo,
                           cp_hi_perm = st_p$cp_upper1,
                           cp_hi_flip = st_f$cp_upper1)
}

tab <- do.call(rbind, rows)
saveRDS(tab, "out/07_size_signflip_summary.rds")
nul <- tab[tab$beta == 0, ]
cat("\n---- verdicts ----\n")
cat(sprintf("  sign-flip null cells with size <= alpha (point estimate)   : %d / %d\n",
            sum(nul$rej_flip <= ALPHA), nrow(nul)))
cat(sprintf("  sign-flip null cells with evidence of over-rejection       : %d / %d\n",
            sum(nul$cp_lo_flip > ALPHA), nrow(nul)))
cat(sprintf("  permutation size rises from homoskedastic to strong het    : %s\n",
            if (nul$rej_perm[1] > nul$rej_perm[3]) "yes" else "NO"))
cat(sprintf("  sign-flip below permutation at the strongest het           : %s\n",
            if (nul$rej_flip[1] <= nul$rej_perm[1]) "yes" else "NO"))
pw <- tab[tab$beta > 0, ]
cat(sprintf("  power at beta = .15, homoskedastic: perm %s%%, flip %s%%  (%s)\n",
            fmt_pct(pw$rej_perm), fmt_pct(pw$rej_flip),
            if (pw$rej_flip < pw$rej_perm)
              "sign-flip is the less powerful, as documented" else
              "UNEXPECTED: sign-flip not less powerful"))
cat("\nfull table:\n")
print(tab, digits = 3, row.names = FALSE)
sink()
