## mwperm_dyadic() -- the two-way (dyadic) design, condition InvA.
##
## Row and column clusters are permuted jointly by a common group order, which
## is exact whenever the error array is separately exchangeable in both
## dimensions (Guo, Toulis & Wang 2026, Assumption 1 and Section 3). This file
## pins what the front end itself owns: the seeded numbers on the shipped data,
## the group order K it derives from the array, how it turns arbitrary cluster
## labels into dense ids, and the one-observation-per-cell contract that makes
## the gather vectors bijections.
##
## The cross-cutting contracts live elsewhere: input validation in
## test-validation.R, invariance and parallel identity in test-equivariance.R,
## the permutation algebra in `lower-level-tests/`. See tests/README.md.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

data(trade_dyadic)

## ---- 1. the seeded anchor on the shipped data -----------------------------
## The documented reference call. tests/golden/baseline.rds is the
## authoritative snapshot; this repeats the headline numbers here so a failure
## names the design that moved rather than "the baseline".
fit <- with(trade_dyadic,
            mwperm_dyadic(y = log_trade, d = log_dist,
                          x = cbind(log_gdp_i, log_gdp_j),
                          row = importer, col = exporter, seed = 1))
stopifnot(inherits(fit, "mwperm"),
          identical(fit$type, "dyadic"),
          identical(fit$pvalue, 0.025),
          abs(fit$estimate - (-0.898503)) < 1e-6,
          abs(fit$conf_int[1] - (-1.251351)) < 1e-5,
          abs(fit$conf_int[2] - (-0.534468)) < 1e-5,
          fit$n_obs == 1600L,
          identical(unname(fit$n_clusters), c(40L, 40L)),
          identical(names(fit$n_clusters), c("row", "col")))
## the true coefficient is -1 (attr "true_coef"); the interval covers it
stopifnot(fit$conf_int[1] <= -1, fit$conf_int[2] >= -1)

## same seed, same answer -- every field, not just the p-value
fit_b <- with(trade_dyadic,
              mwperm_dyadic(y = log_trade, d = log_dist,
                            x = cbind(log_gdp_i, log_gdp_j),
                            row = importer, col = exporter, seed = 1))
stopifnot(isTRUE(same_fit(fit, fit_b)))

## the placebo column carries no signal: the test must not reject it
plc <- with(trade_dyadic,
            mwperm_dyadic(y = log_trade, d = placebo,
                          x = cbind(log_gdp_i, log_gdp_j),
                          row = importer, col = exporter, seed = 1,
                          conf_int = FALSE))
stopifnot(identical(plc$pvalue, 0.975))

## ---- 2. the group order K comes from the permuted dimensions --------------
## Default K = min(n_row, n_col) - 1, capped at 199, and n_perm = K + 1.
stopifnot(fit$K == 39L, fit$n_perm == 40L,
          identical(fit$resolution, 1 / 40), identical(fit$p_floor, 1 / 40),
          length(fit$pvalues_rep) == fit$n_reps, fit$n_reps == 10L)

set.seed(11)
g69 <- expand.grid(i = 1:6, j = 1:9)              # non-square: K follows 6
y69 <- rnorm(nrow(g69))
d69 <- rnorm(nrow(g69))
f69 <- mwperm_dyadic(y69, d69, row = g69$i, col = g69$j, seed = 1,
                     conf_int = FALSE)
stopifnot(f69$K == 5L, identical(unname(f69$n_clusters), c(6L, 9L)))

## an explicit K is honoured, and K + 1 may not exceed the smaller dimension
f69k <- mwperm_dyadic(y69, d69, row = g69$i, col = g69$j, K = 3L, seed = 1,
                      conf_int = FALSE)
stopifnot(f69k$K == 3L, f69k$n_perm == 4L)
expect_err(mwperm_dyadic(y69, d69, row = g69$i, col = g69$j, K = 6L), "K + 1")

## ---- 3. cluster labels of any type give the same dense coding -------------
## .dense_id() maps labels to 1..n_cluster; the fit must not depend on whether
## the user supplied integers, characters, factors or gappy codes.
set.seed(101)
g6 <- expand.grid(i = 1:6, j = 1:6)
N6 <- nrow(g6)
d6 <- rnorm(6)[g6$i] + rnorm(N6)
y6 <- rnorm(6)[g6$i] + rnorm(6)[g6$j] + 0.4 * d6 + rnorm(N6)

f_int <- mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 1,
                       conf_int = FALSE)
f_chr <- mwperm_dyadic(y6, d6, row = letters[g6$i], col = g6$j, seed = 1,
                       conf_int = FALSE)
f_gap <- mwperm_dyadic(y6, d6, row = c(10L, 20L, 35L, 40L, 70L, 99L)[g6$i],
                       col = g6$j, seed = 1, conf_int = FALSE)
stopifnot(identical(f_chr$pvalue, f_int$pvalue),
          identical(f_gap$pvalue, f_int$pvalue))

## unused factor levels are dropped, not counted as empty clusters
f_fac <- mwperm_dyadic(y6, d6, row = factor(g6$i, levels = 1:10), col = g6$j,
                       seed = 1, conf_int = FALSE)
stopifnot(f_fac$n_clusters[["row"]] == 6L,
          identical(f_fac$pvalue, f_int$pvalue))

## ---- 4. exactly one observation per (row, col) cell -----------------------
## A repeated cell would make the gather vector many-to-one: match() resolves
## the duplicate to its first hit, so the statistic would be computed on
## duplicated rows with no NA and no warning. It must be an error instead.
expect_err(mwperm_dyadic(c(y6, 1), c(d6, 1), row = c(g6$i, 1),
                         col = c(g6$j, 1)),
           "one observation per (row, col) cell")
## replicated cells are a different design, and the message says so
stopifnot(grepl("mwperm_layout",
                msg_of(mwperm_dyadic(c(y6, y6), c(d6, d6), row = c(g6$i, g6$i),
                                     col = c(g6$j, g6$j)))))

## An INCOMPLETE array is refused up front, from the cell count, before any
## permutation is drawn. Catching it inside the gather-vector builder (when a
## draw happens to reach an unobserved cell) is in principle draw-dependent: a
## group that mapped the observed set onto itself would slip through and run on
## a non-rectangular design with K set from the full id range.
g_inc <- g6[g6$i != g6$j, ]                        # diagonal deleted
for (sd in 1:5) {
  m_inc <- msg_of(mwperm_dyadic(y6[g6$i != g6$j], d6[g6$i != g6$j],
                                row = g_inc$i, col = g_inc$j, seed = sd,
                                n_reps = 1, conf_int = FALSE))
  stopifnot(!is.na(m_inc), grepl("complete", m_inc, fixed = TRUE),
            grepl("expected 36", m_inc, fixed = TRUE),
            grepl("found 30", m_inc, fixed = TRUE),
            grepl("mwperm_missing()", m_inc, fixed = TRUE),
            !grepl("Permutation maps", m_inc, fixed = TRUE))
}

## ---- 5. the coefficient name follows the supplied column ------------------
stopifnot(identical(fit$d_names, "log_dist"))
f_nm <- mwperm_dyadic(y6, cbind(treat = d6), row = g6$i, col = g6$j, seed = 1,
                      conf_int = FALSE)
stopifnot(identical(f_nm$d_names, "treat"))
f_ex <- mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 1,
                      conf_int = FALSE)
stopifnot(identical(f_ex$d_names, "d6"))          # deparse(substitute(d))

passed("test-dyadic.R")
