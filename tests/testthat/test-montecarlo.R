## Finite-sample validity by Monte Carlo (slow; MWPERM_SLOW_TESTS=true).
##
## Theorem 1 of Guo, Toulis & Wang (2026): under H0, P(pval <= alpha | X, D)
## <= alpha for EVERY alpha. With a single repetition (the configuration the
## theorem covers) the p-value lives on {1, ..., K + 1}/(K + 1), so the
## claim is checked at every atom: the empirical CDF at a = j/(K + 1) must
## not exceed a + 3 binomial standard errors (expect_valid_ecdf()). Designs:
## two-way clustered heavy-tailed errors (t_3, Cauchy), the node-level
## Moulton regressor of the JSS draft's Section 5, an MCAR mask for
## Procedure 2 (Theorem 4, conditional on the mask), and a panel with a
## random-walk common trend and AR(1) errors (InvB). Then the draft's
## Section 3.2 negative control (the three-way test over-rejects on a
## trending panel; the panel test does not) and a power check. Every
## simulation is seeded; each design uses K + 1 >= 20 so that rejection at
## 0.05 is attainable and a small rejection rate means something.

n_sim <- 1000L

sim_p <- function(n_sim, seed0, one) vapply(seq_len(n_sim), function(s)
  one(seed0 + s), numeric(1))

test_that("dyadic, two-way random effects, t(3) errors: valid at every atom", {
  skip_if_not_slow()
  p <- sim_p(n_sim, 10000, function(s) {
    dy <- make_dyadic(20, 20, beta = 0, seed = s,
                      err = function(k) stats::rt(k, 3))
    mwperm_dyadic(dy$y, dy$d, x = dy$x, row = dy$i, col = dy$j, n_reps = 1,
                  seed = s, conf_int = FALSE)$pvalue
  })
  mc_record("dyadic_t3", expect_valid_ecdf(p, 19))
})

test_that("dyadic, two-way random effects, Cauchy errors: valid at every atom", {
  skip_if_not_slow()
  p <- sim_p(n_sim, 20000, function(s) {
    dy <- make_dyadic(20, 20, beta = 0, seed = s,
                      err = function(k) stats::rcauchy(k))
    mwperm_dyadic(dy$y, dy$d, x = dy$x, row = dy$i, col = dy$j, n_reps = 1,
                  seed = s, conf_int = FALSE)$pvalue
  })
  mc_record("dyadic_cauchy", expect_valid_ecdf(p, 19))
})

test_that("the Moulton design (draft Sec. 5.1): IPT valid, naive OLS is not", {
  skip_if_not_slow()
  ## n = 25, phi = 0.3, eps = sqrt(phi)(eta_i + xi_j) + sqrt(1 - phi) u_ij
  ## with t_4 / sqrt(2) components, and a node-level d_ij = d_i.
  one <- function(s) with_seed(s, {
    n <- 25
    phi <- 0.3
    g <- expand.grid(i = seq_len(n), j = seq_len(n))
    rt4 <- function(k) stats::rt(k, df = 4) / sqrt(2)
    eta <- rt4(n)
    xi <- rt4(n)
    u <- rt4(n^2)
    d <- stats::rnorm(n)[g$i]
    y <- sqrt(phi) * (eta[g$i] + xi[g$j]) + sqrt(1 - phi) * u
    c(ipt = mwperm_dyadic(y, d, row = g$i, col = g$j, n_reps = 1, seed = s,
                          conf_int = FALSE)$pvalue,
      ols = summary(stats::lm(y ~ d))$coefficients["d", 4])
  })
  res <- vapply(30000 + seq_len(n_sim), one, numeric(2))
  mc_record("moulton", expect_valid_ecdf(res["ipt", ], 24))
  expect_gt(mean(res["ols", ] <= 0.05), 0.2)       # the naive test fails
})

test_that("missing data, MCAR mask (Procedure 2, Theorem 4): valid", {
  skip_if_not_slow()
  ## One MCAR mask on a 40 x 40 array, held fixed (Theorem 4 conditions on
  ## M); its greedy biclique is 27 x 25, so K = 19 is attainable.
  keep <- with_seed(7, stats::runif(1600) < 0.96)
  p <- sim_p(n_sim, 40000, function(s) {
    dy <- make_dyadic(40, 40, beta = 0, seed = s,
                      err = function(k) stats::rt(k, 3))[keep, ]
    mwperm_missing(dy$y, dy$d, x = dy$x, row = dy$i, col = dy$j,
                   min_block = 20, K = 19, n_reps = 1, seed = s,
                   conf_int = FALSE)$pvalue
  })
  mc_record("missing_mcar", expect_valid_ecdf(p, 19))
})

test_that("panel with a random-walk trend and AR(1) errors (InvB): valid", {
  skip_if_not_slow()
  p <- sim_p(n_sim, 50000, function(s) {
    pn <- make_panel(20, 20, T = 4, beta = 0, seed = s, trend_sd = 3,
                     rho = 0.7)
    mwperm_panel(pn$y, pn$d, x = pn$x, row = pn$i, col = pn$j, time = pn$t,
                 n_reps = 1, seed = s, conf_int = FALSE)$pvalue
  })
  mc_record("panel_trend", expect_valid_ecdf(p, 19))
})

test_that("negative control (draft Sec. 3.2): three-way over-rejects, panel holds", {
  skip_if_not_slow()
  ## 20 x 20 pairs over 20 periods, staggered treatment d = 1{t >= start_ij},
  ## a random-walk common time effect, beta = 0. InvA is false (time is not
  ## exchangeable); InvB holds.
  one <- function(s) with_seed(s, {
    n <- 20
    Tt <- 20
    g <- expand.grid(i = seq_len(n), j = seq_len(n), t = seq_len(Tt))
    start <- matrix(sample.int(Tt + 1L, n * n, replace = TRUE), n, n)
    d <- as.numeric(g$t >= start[cbind(g$i, g$j)])
    zeta <- cumsum(stats::rnorm(Tt))
    y <- stats::rnorm(n)[g$i] + stats::rnorm(n)[g$j] + zeta[g$t] +
      stats::rnorm(nrow(g))
    c(threeway = mwperm_threeway(y, d, id1 = g$i, id2 = g$j, id3 = g$t,
                                 n_reps = 1, seed = s,
                                 conf_int = FALSE)$pvalue,
      panel = mwperm_panel(y, d, row = g$i, col = g$j, time = g$t,
                           n_reps = 1, seed = s, conf_int = FALSE)$pvalue)
  })
  res <- vapply(60000 + seq_len(200), one, numeric(2))
  rej <- rowMeans(res <= 0.05)
  mc_record("negative_control", data.frame(
    atom = 0.05, ecdf = rej, bound = 0.05 + 3 * sqrt(0.05 * 0.95 / 200),
    test = names(rej)))
  expect_gt(rej[["threeway"]], 0.30)                 # draft: 0.629
  expect_lte(rej[["panel"]], 0.05 + 3 * sqrt(0.05 * 0.95 / 200))
})

test_that("power: a large beta is rejected almost always", {
  skip_if_not_slow()
  p <- sim_p(200, 70000, function(s) {
    dy <- make_dyadic(20, 20, beta = 0.6, seed = s)
    mwperm_dyadic(dy$y, dy$d, x = dy$x, row = dy$i, col = dy$j, n_reps = 1,
                  seed = s, conf_int = FALSE)$pvalue
  })
  mc_record("power_beta_0.6", data.frame(atom = 0.05,
                                         ecdf = mean(p <= 0.05), bound = NA))
  expect_gte(mean(p <= 0.05), 0.95)
})
