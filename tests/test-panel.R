## mwperm_panel() -- the panel design, condition InvB.
##
## The SAME row/column relabelling (pi, sigma) is applied in every period and
## the time index is held fixed, which is exact under within-period
## exchangeability even when an ARBITRARY common time trend zeta_t is present
## (Guo, Toulis & Wang 2026, Section 6.2). That trend robustness is the reason
## the design exists, so section 2 asserts it directly rather than through a
## simulation. `time_fe` adds period dummies -- invariant to the within-period
## permutation, hence valid -- to de-bias the point estimate.
##
## Cross-cutting contracts live elsewhere; see tests/README.md.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

data(trade_panel)

## ---- 1. the seeded anchor on the shipped data -----------------------------
fit <- with(trade_panel,
            mwperm_panel(y = log_trade, d = fta,
                         x = cbind(log_gdp_i, log_gdp_j),
                         row = importer, col = exporter, time = year,
                         seed = 1))
stopifnot(inherits(fit, "mwperm"),
          identical(fit$type, "panel"),
          abs(fit$pvalue - 0.0454545454545455) < 1e-12,
          abs(fit$estimate - 0.677407) < 1e-6,
          abs(fit$conf_int[1] - 0.442022) < 1e-5,
          abs(fit$conf_int[2] - 0.880272) < 1e-5,
          identical(unname(fit$n_clusters), c(22L, 22L, 6L)),
          identical(names(fit$n_clusters), c("row", "col", "time")))
fit_b <- with(trade_panel,
              mwperm_panel(y = log_trade, d = fta,
                           x = cbind(log_gdp_i, log_gdp_j),
                           row = importer, col = exporter, time = year,
                           seed = 1))
stopifnot(isTRUE(same_fit(fit, fit_b)))

## ---- 2. InvB: an arbitrary common time trend changes nothing --------------
## With `time_fe = TRUE` the period dummies span every zeta_t, so adding one to
## the outcome shifts y inside col(X): by FWL the residualized problem is
## untouched. The p-value must be BIT-identical, the estimate and interval
## equal to machine epsilon. No assumption is made about the trend's shape --
## this one is neither smooth nor monotone.
tt <- as.integer(factor(trade_panel$year))
zeta <- c(0, 3.1, -2.4, 7.0, 1.5, -4.2)[tt]
fit_z <- with(trade_panel,
              mwperm_panel(y = log_trade + zeta, d = fta,
                           x = cbind(log_gdp_i, log_gdp_j),
                           row = importer, col = exporter, time = year,
                           seed = 1))
stopifnot(identical(fit$pvalue, fit_z$pvalue),
          identical(fit$pvalues_rep, fit_z$pvalues_rep),
          max(abs(fit$estimate - fit_z$estimate)) < 1e-10,
          max(abs(fit$conf_int - fit_z$conf_int)) < 1e-10)

## Without the dummies the trend is not projected out, so the estimate moves --
## which is exactly why `time_fe = TRUE` is the default. On these data the true
## coefficient is 0.5 (attr "true_coef"); the de-biased estimate is closer.
no_fe <- with(trade_panel,
              mwperm_panel(y = log_trade, d = fta,
                           x = cbind(log_gdp_i, log_gdp_j),
                           row = importer, col = exporter, time = year,
                           seed = 1, time_fe = FALSE, conf_int = FALSE))
stopifnot(abs(no_fe$estimate - 0.961030) < 1e-6,
          abs(fit$estimate - 0.5) < abs(no_fe$estimate - 0.5))

## ---- 3. K is set by the permuted dimensions, never by time ----------------
## Time is held fixed, so the number of periods must not enter the cap: a
## 6 x 6 x 4 panel has K = 5, not 3.
set.seed(21)
gp <- expand.grid(i = 1:6, j = 1:6, t = 1:4)
Np <- nrow(gp)
dp <- rnorm(Np)
yp <- rnorm(6)[gp$i] + rnorm(6)[gp$j] + cumsum(rnorm(4))[gp$t] + 0.3 * dp +
  rnorm(Np)
fp <- mwperm_panel(yp, dp, row = gp$i, col = gp$j, time = gp$t, seed = 3,
                   conf_int = FALSE)
stopifnot(fp$K == 5L, fp$n_perm == 6L)
expect_err(mwperm_panel(yp, dp, row = gp$i, col = gp$j, time = gp$t, K = 6L),
           "K + 1")

## ---- 4. row-order invariance ----------------------------------------------
## Observations are gathered by (row, col, time) cell, so the order in which
## the user happens to supply them cannot matter.
o <- sample(Np)
fp_o <- mwperm_panel(yp[o], dp[o], row = gp$i[o], col = gp$j[o],
                     time = gp$t[o], seed = 3, conf_int = FALSE)
stopifnot(identical(fp$pvalue, fp_o$pvalue),
          abs(fp$estimate - fp_o$estimate) < 1e-10)

## ---- 5. the array must be complete and balanced ---------------------------
## InvB permutes the same (pi, sigma) in every period, which is only defined
## when every (row, col) pair is observed in every period.
expect_err(mwperm_panel(yp[-1], dp[-1], row = gp$i[-1], col = gp$j[-1],
                        time = gp$t[-1]),
           "complete")
expect_err(mwperm_panel(c(yp, 1), c(dp, 1), row = c(gp$i, 1),
                        col = c(gp$j, 1), time = c(gp$t, 1)),
           "at most once")

passed("test-panel.R")
