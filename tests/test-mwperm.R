## Tests for the user-facing tests: dyadic validity, CI behaviour,
## biclique finder, and the missing-data path. Kept small so they run fast
## under R CMD check. Plain base-R (no testthat).

library(mwperm)

## ---- a small reproducible dyadic data set --------------------------------
make_dyadic <- function(n = 16, beta = 1, seed = 1) {
  set.seed(seed)
  g <- expand.grid(i = seq_len(n), j = seq_len(n))
  g <- g[g$i != g$j, ]                       # drop diagonal -> some missing
  eta <- rnorm(n); xi <- rnorm(n)
  g$d <- rnorm(nrow(g))
  g$z <- rnorm(nrow(g))                      # an irrelevant covariate
  g$y <- beta * g$d + 0.5 * g$z + eta[g$i] + xi[g$j] +
    0.5 * rt(nrow(g), df = 4)
  g
}

## ---- (1) dyadic test on a complete array ---------------------------------
## Need >= 20 clusters per dimension so that 1/(K+1) <= 0.05 and a 95% CI
## is attainable; use 22x22.
g <- expand.grid(i = 1:22, j = 1:22)         # complete 22x22
set.seed(7)
eta <- rnorm(22); xi <- rnorm(22)
g$d <- rnorm(nrow(g)); g$z <- rnorm(nrow(g))
g$y <- 1.0 * g$d + 0.5 * g$z + eta[g$i] + xi[g$j] + 0.5 * rt(nrow(g), 4)

fit <- with(g, mwperm_dyadic(y = y, d = d, x = z, row = i, col = j,
                             n_reps = 5, seed = 11))
stopifnot(inherits(fit, "mwperm"))
stopifnot(fit$pvalue >= 0, fit$pvalue <= 1)
stopifnot(abs(fit$estimate - 1.0) < 0.5)      # estimate near the truth
stopifnot(fit$pvalue < 0.10)                  # genuine signal detected

## CI must bracket the point estimate and exclude 0 here
stopifnot(!is.null(fit$conf_int))
stopifnot(fit$conf_int[1] <= fit$estimate, fit$estimate <= fit$conf_int[2])
stopifnot(fit$conf_int[1] > 0)

## ---- (2) placebo (true coefficient 0) should rarely reject ---------------
fitp <- with(g, mwperm_dyadic(y = y, d = z, x = d, row = i, col = j,
                              conf_int = TRUE, n_reps = 5, seed = 13))
## z has a real effect too here, so instead test a pure-noise regressor:
set.seed(99); g$noise <- rnorm(nrow(g))
fitn <- with(g, mwperm_dyadic(y = y, d = noise, x = cbind(d, z),
                              row = i, col = j, n_reps = 5, seed = 17))
stopifnot(fitn$pvalue > 0.10)                 # no false signal
stopifnot(fitn$conf_int[1] < 0, fitn$conf_int[2] > 0)  # CI covers 0

## ---- (3) summary()/confint() contract ------------------------------------
s <- summary(fit)
stopifnot(is.data.frame(s))
stopifnot(all(c("term", "estimate", "conf_low", "conf_high", "p_value") %in%
              names(s)))
ci <- confint(fit)
stopifnot(is.matrix(ci), nrow(ci) == 1L)
stopifnot(abs(ci[1, 1] - fit$conf_int[1]) < 1e-8)

## confint() with a non-stored level must error (cannot be re-derived)
err <- tryCatch({ confint(fit, level = 0.90); FALSE }, error = function(e) TRUE)
stopifnot(err)

## ---- (4) coarse resolution -> no CI when 1/(K+1) > alpha -----------------
small <- expand.grid(i = 1:6, j = 1:6)        # K defaults to 5 -> min p = 1/6
set.seed(3)
ei <- rnorm(6); ej <- rnorm(6)
small$d <- rnorm(36)
small$y <- small$d + ei[small$i] + ej[small$j] + rnorm(36)
fits <- with(small, mwperm_dyadic(y = y, d = d, row = i, col = j,
                                  alpha = 0.05, seed = 1))
stopifnot(is.null(fits$conf_int))             # 1/6 = 0.167 > 0.05
stopifnot(!is.null(fits$note))

## ---- (5) biclique finder: all-ones and disjoint --------------------------
gm <- make_dyadic(n = 16, seed = 5)
blocks <- find_bicliques(gm$i, gm$j, min_block = 3)
stopifnot(length(blocks) >= 1)
obs <- paste(gm$i, gm$j, sep = "-")
seen_rows <- integer(0); seen_cols <- integer(0)
for (b in blocks) {
  ## every cell in rows x cols must be observed (all-ones)
  cells <- expand.grid(r = b$rows, c = b$cols)
  stopifnot(all(paste(cells$r, cells$c, sep = "-") %in% obs))
  ## blocks disjoint in rows and in cols
  stopifnot(length(intersect(seen_rows, b$rows)) == 0)
  stopifnot(length(intersect(seen_cols, b$cols)) == 0)
  seen_rows <- c(seen_rows, b$rows); seen_cols <- c(seen_cols, b$cols)
}

## ---- (6) missing-data path returns a valid result ------------------------
fitm <- with(gm, mwperm_missing(y = y, d = d, x = z, row = i, col = j,
                                min_block = 3, n_reps = 3, seed = 9))
stopifnot(inherits(fitm, "mwperm"))
stopifnot(fitm$pvalue >= 0, fitm$pvalue <= 1)
stopifnot(fitm$n_obs <= nrow(gm))             # only fully-observed cells used
stopifnot(grepl("missing", fitm$type, fixed = TRUE))

## ---- (7) non-zero null: testing H0 at the truth must not reject -----------
## exercises the permuted-D slope used in the b-statistic when beta_null != 0
## (with conf_int = FALSE, the fast path that skips W only when the null is 0).
fit_true <- with(g, mwperm_dyadic(y = y, d = d, x = z, row = i, col = j,
                                  beta_null = 1.0, conf_int = FALSE,
                                  n_reps = 5, seed = 11))
fit_zero <- with(g, mwperm_dyadic(y = y, d = d, x = z, row = i, col = j,
                                  beta_null = 0.0, conf_int = FALSE,
                                  n_reps = 5, seed = 11))
stopifnot(fit_true$pvalue > fit_zero$pvalue)   # truth (1.0) much harder to reject than 0
stopifnot(fit_true$pvalue > 0.10)              # do not reject the true value

## ---- (8) layout L0 balancing (Section 6.3) -------------------------------
set.seed(31)
lr <- integer(0); lc <- integer(0)
for (i in 1:5) for (j in 1:5) {
  m_ij <- sample(3:7, 1)
  lr <- c(lr, rep(i, m_ij)); lc <- c(lc, rep(j, m_ij))
}
ld <- rnorm(length(lr))
ly <- 0.5 * ld + rnorm(5)[lr] + rnorm(5)[lc] + 0.3 * rnorm(length(lr))
fitL <- mwperm_layout(y = ly, d = ld, row = lr, col = lc, L0 = 4L,
                      n_reps = 3, seed = 5, conf_int = FALSE)
n_kept <- sum(table(interaction(lr, lc, drop = TRUE)) >= 4L)
stopifnot(fitL$n_obs == 4L * n_kept)           # every kept cell downsampled to L0
stopifnot(fitL$K == 3L)                        # K defaults to L0 - 1
stopifnot(grepl("balanced", fitL$type, fixed = TRUE))
## reproducible given the seed
fitL2 <- mwperm_layout(y = ly, d = ld, row = lr, col = lc, L0 = 4L,
                       n_reps = 3, seed = 5, conf_int = FALSE)
stopifnot(identical(fitL$pvalue, fitL2$pvalue), identical(fitL$n_obs, fitL2$n_obs))

## ---- (9) exact biclique method: valid, disjoint, >= greedy ---------------
## small structured mask where the exact search completes within budget
set.seed(51)
em <- expand.grid(i = 1:10, j = 1:10)
em <- em[runif(nrow(em)) < 0.8, ]            # 80% observed
bex <- find_bicliques(em$i, em$j, min_block = 2, method = "exact")
bgr <- find_bicliques(em$i, em$j, min_block = 2, method = "greedy")
obs9 <- paste(em$i, em$j, sep = "-")
area <- function(bs) sum(vapply(bs, function(b) length(b$rows) * length(b$cols), 0))
sr <- integer(0); sc <- integer(0)
for (b in bex) {
  cells <- expand.grid(r = b$rows, c = b$cols)
  stopifnot(all(paste(cells$r, cells$c, sep = "-") %in% obs9))   # all-ones
  stopifnot(length(intersect(sr, b$rows)) == 0,
            length(intersect(sc, b$cols)) == 0)                  # disjoint
  sr <- c(sr, b$rows); sc <- c(sc, b$cols)
}
stopifnot(area(bex) >= area(bgr))            # exact never worse than greedy

## mwperm_missing with block_method = "exact" returns a valid result
fme <- with(gm, mwperm_missing(y = y, d = d, x = z, row = i, col = j,
                               min_block = 3, n_reps = 3, seed = 9,
                               block_method = "exact"))
stopifnot(inherits(fme, "mwperm"))
stopifnot(fme$pvalue >= 0, fme$pvalue <= 1)

## ---- (10) multi-coefficient joint confidence region (d > 1) --------------
fit2 <- with(g, mwperm_dyadic(y = y, d = cbind(d, z), row = i, col = j,
                              conf_int = TRUE, n_reps = 3, seed = 11))
stopifnot(is.null(fit2$conf_int))                 # no scalar interval for d > 1
stopifnot(!is.null(fit2$conf_box), is.matrix(fit2$conf_box))
stopifnot(dim(fit2$conf_box)[1] == 2L, dim(fit2$conf_box)[2] == 2L)
stopifnot(nrow(fit2$conf_region) >= 1L)
## point estimates lie within the marginal box
stopifnot(all(fit2$estimate >= fit2$conf_box[1, ] &
              fit2$estimate <= fit2$conf_box[2, ]))
## confint() returns one row per coefficient; summary() fills both rows
ci2 <- confint(fit2)
stopifnot(is.matrix(ci2), nrow(ci2) == 2L, ncol(ci2) == 2L)
s2 <- summary(fit2)
stopifnot(nrow(s2) == 2L, all(!is.na(s2$conf_low)), all(!is.na(s2$conf_high)))
## a vector null is accepted and tested jointly
fitv <- with(g, mwperm_dyadic(y = y, d = cbind(d, z), row = i, col = j,
                              beta_null = c(1.0, 0.5), conf_int = FALSE, seed = 11))
stopifnot(fitv$pvalue >= 0, fitv$pvalue <= 1)
## coarse resolution (d > 1): no region, explanatory note
fit2s <- with(small, mwperm_dyadic(y = y, d = cbind(d, d^2), row = i, col = j, seed = 1))
stopifnot(is.null(fit2s$conf_box), length(fit2s$note) > 0)

cat("test-mwperm.R: all checks passed\n")
