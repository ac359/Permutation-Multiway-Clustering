## Exact invariances of Procedure 1 (Guo, Toulis & Wang 2026) and the
## equivariance of its inverted confidence set.
##
## For a fixed seed the permutation groups are fixed, so each invariance is a
## deterministic identity between two fits: the whole `pvalues_rep` vector is
## compared (tolerance 1e-10, i.e. no indicator in Eq. 10 may flip), and the
## confidence sets to 1e-8 relative.
##   * V_k' X = 0, so y -> y + X gamma changes no statistic (Procedure 1,
##     step 1), however large gamma is;
##   * a_k and b_k are both scaled by |c| under y -> c y, and by |c| under
##     D -> c D, so their comparison in Eq. (10) is unchanged;
##   * the test of beta = b IS the test of beta = 0 on y - D b (step 3);
##   * with period dummies in X (time_fe = TRUE), a common trend zeta_t lies in
##     col(X): condition InvB of GTW Section 6.2.
## Complements tests/test-equivariance.R, which checks row order and small
## nuisance shifts on a different fixture.

fx <- make_dyadic(20, 20, beta = 0.3, seed = 31)
fit_of <- function(y, d = fx$d, x = fx$x, ...)
  mwperm_dyadic(y, d, x = x, row = fx$i, col = fx$j, n_reps = 3, seed = 41,
                ...)
base <- fit_of(fx$y)

test_that("the fixture supports a 95% set and has signal", {
  expect_identical(base$K, 19L)
  expect_false(is.null(base$conf_int))
})

test_that("y + X gamma leaves every p-value and the set unchanged", {
  X <- cbind(1, fx$x)
  for (gamma in list(c(1e3, -2e3), c(-5e4, 7.5e3))) {
    f <- fit_of(fx$y + drop(X %*% gamma))
    expect_equal(f$pvalues_rep, base$pvalues_rep, tolerance = 1e-10)
    expect_equal(f$conf_int, base$conf_int, tolerance = 1e-8)
    expect_equal(f$conf_set, base$conf_set, tolerance = 1e-8)
  }
})

test_that("rescaling y: p-values invariant, CI(c y) = c CI(y)", {
  for (c in c(0.001, 7.5, 1e4)) {
    f <- fit_of(c * fx$y)
    expect_equal(f$pvalues_rep, base$pvalues_rep, tolerance = 1e-10)
    expect_equal(f$conf_int, c * base$conf_int, tolerance = 1e-8)
  }
  f <- fit_of(-2 * fx$y)                  # a negative scale reverses the set
  expect_equal(f$pvalues_rep, base$pvalues_rep, tolerance = 1e-10)
  expect_equal(f$conf_int, rev(-2 * base$conf_int), tolerance = 1e-8)
})

test_that("rescaling D: p-values invariant, CI(y; c D) = CI(y; D) / c", {
  for (c in c(0.01, 3, -4)) {
    f <- fit_of(fx$y, d = c * fx$d)
    expect_equal(f$pvalues_rep, base$pvalues_rep, tolerance = 1e-10)
    expect_equal(sort(f$conf_int), sort(base$conf_int / c), tolerance = 1e-8)
  }
})

test_that("testing beta = b on y is testing beta = 0 on y - D b", {
  for (b in c(-0.4, 0.25, 0.3, 2)) {
    f1 <- fit_of(fx$y, beta_null = b, conf_int = FALSE)
    f0 <- fit_of(fx$y - b * fx$d, conf_int = FALSE)
    expect_equal(f1$pvalues_rep, f0$pvalues_rep, tolerance = 1e-10)
  }
})

test_that("CI(y + c D) = CI(y) + c", {
  for (c in c(-3, 0.5, 40)) {
    f <- fit_of(fx$y + c * fx$d)
    expect_equal(f$conf_int, base$conf_int + c, tolerance = 1e-8)
    expect_equal(f$conf_set, base$conf_set + c, tolerance = 1e-8)
  }
})

test_that("panel, time_fe = TRUE: an arbitrary common trend zeta_t is invisible", {
  ## d is correlated with t, so without period dummies the trend WOULD move
  ## the test; with them it cannot.
  pn <- make_panel(20, 20, T = 3, beta = 0.2, seed = 32)
  pfit <- function(y, time_fe = TRUE)
    mwperm_panel(y, pn$d, x = pn$x, row = pn$i, col = pn$j, time = pn$t,
                 n_reps = 3, seed = 42, time_fe = time_fe)
  b0 <- pfit(pn$y)
  zeta <- c(250, -1e3, 3.3e4)                        # arbitrary, huge
  b1 <- pfit(pn$y + zeta[pn$t])
  expect_equal(b1$pvalues_rep, b0$pvalues_rep, tolerance = 1e-10)
  expect_equal(b1$conf_int, b0$conf_int, tolerance = 1e-8)
  ## the same trend does change the fit when it is NOT partialled out
  n0 <- pfit(pn$y, time_fe = FALSE)
  n1 <- pfit(pn$y + zeta[pn$t], time_fe = FALSE)
  expect_false(isTRUE(all.equal(n0$conf_int, n1$conf_int)))
})
