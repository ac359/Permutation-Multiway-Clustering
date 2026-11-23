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
##     col(X): condition InvB of GTW Section 6.2;
##   * the projector depends on col([X | X_k]) alone, so the units and the
##     origin of a nuisance column cannot move a statistic either.
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

test_that("units and origin of a nuisance column: p-values and set unchanged", {
  ## col([X | X_k]) does not change when a column of x is rescaled, or
  ## shifted (the intercept is in X). The p-values are compared EXACTLY (each
  ## is a count of comparisons) and the set to a tolerance. The second
  ## nuisance column w carries signal, so a w that is silently not
  ## partialled out shows: on the dyadic fixture, whose p-value sits at its
  ## floor 1/(K+1), in conf_int; in the sign-flip and panel fits below
  ## also in pvalues_rep.
  w <- with_seed(37, exp(stats::rnorm(nrow(fx))))
  y <- fx$y + 0.5 * w
  x0 <- cbind(fx$x, w)
  b0 <- fit_of(y, x = x0)
  for (x1 in list(cbind(fx$x, w * 1e9), cbind(fx$x, w * 1e-9),
                  cbind(fx$x, w + 1e6), cbind(fx$x * 1e7, w))) {
    f <- fit_of(y, x = x1)
    expect_identical(f$pvalues_rep, b0$pvalues_rep)
    expect_equal(f$conf_int, b0$conf_int, tolerance = 1e-6)
  }
  ## the sign-flip group takes the same projection
  het <- function(x) mwperm_dyadic_het(y, fx$d, x = x, row = fx$i,
                                       col = fx$j, n_reps = 3, seed = 41)
  h0 <- het(x0)
  h1 <- het(cbind(fx$x, w * 1e9))
  expect_identical(h1$pvalues_rep, h0$pvalues_rep)
  expect_equal(h1$conf_int, h0$conf_int, tolerance = 1e-6)
  ## and the panel design with period dummies (time_fe = TRUE)
  pn <- make_panel(20, 20, T = 3, beta = 0.2, seed = 32)
  pw <- with_seed(38, exp(stats::rnorm(nrow(pn))))
  pfit <- function(x) mwperm_panel(pn$y + pw, pn$d, x = x, row = pn$i,
                                   col = pn$j, time = pn$t, n_reps = 3,
                                   seed = 42)
  p0 <- pfit(cbind(pn$x, pw))
  p1 <- pfit(cbind(pn$x, pw * 1e9))
  expect_identical(p1$pvalues_rep, p0$pvalues_rep)
  expect_equal(p1$conf_int, p0$conf_int, tolerance = 1e-6)
})

test_that("an element that fixes d and x keeps p = 1 on the QR route", {
  ## A layout whose d and x are constant within cells: every element fixes
  ## D and X, so a_k == b_k exactly and p = 1 whichever route the
  ## projection takes (the x below is in large units, which can send it to
  ## the QR, whose rows for identical inputs may differ in the last bit).
  dat <- with_seed(101, {
    n <- 6L
    L <- 21L
    row <- rep(rep(seq_len(n), each = n), each = L)
    col <- rep(rep(seq_len(n), times = n), each = L)
    cid <- (row - 1L) * n + col
    dc <- stats::rnorm(n * n)
    xc <- stats::rnorm(n * n)
    list(y = 0.3 * dc[cid] + 0.5 * xc[cid] + stats::rnorm(length(cid)),
         d = dc[cid], x = xc[cid] * 1e5, row = row, col = col)
  })
  f <- suppressWarnings(mwperm_layout(dat$y, dat$d, dat$x, row = dat$row,
                                      col = dat$col, seed = 1))
  expect_identical(f$pvalue, 1)
  expect_identical(f$conf_int, c(-Inf, Inf))
})

test_that("the origin of the covariate of interest changes nothing", {
  ## The intercept is in X, so the test of d + c is the test of d; the
  ## engine centres a d far from zero, where its origin would otherwise
  ## cost accuracy (and, from about 1e8 sd, called it constant).
  for (k in c(1e5, 1e8)) {
    f <- fit_of(fx$y, d = fx$d + k * stats::sd(fx$d))
    expect_identical(f$pvalues_rep, base$pvalues_rep)
    expect_equal(f$conf_int, base$conf_int, tolerance = 1e-6)
  }
  ## a d whose only variation beyond x is the rounding of its own origin
  ## still has none: a warning and p = 1
  x1 <- as.matrix(fx$x)[, 1]
  w <- character(0)
  f0 <- withCallingHandlers(
    fit_of(fx$y, d = x1 + 1e9 * stats::sd(x1), conf_int = FALSE),
    warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    })
  expect_true(any(grepl("no variation", w, fixed = TRUE)))
  expect_identical(f0$pvalue, 1)
})

test_that("adding a constant to y changes nothing, even a very large one", {
  ## The intercept absorbs the constant in exact arithmetic. In floating point
  ## an outcome far from zero multiplies the rounding left in the residualized
  ## regressor by its origin, so the engine centres an outcome whose mean
  ## exceeds 1e3 standard deviations (y + 1e9 used to move the interval by
  ## 8e-5 relative, with no warning). What the caller's own arithmetic rounded
  ## away cannot come back, hence the 1e-7 tolerance at 1e9.
  dy <- make_dyadic(20, 20, seed = 61)
  f <- function(y) mwperm_dyadic(y, dy$d, x = dy$x, row = dy$i, col = dy$j,
                                 seed = 3)
  b <- f(dy$y)
  for (h in c(1e6, 1e9)) {
    g <- f(dy$y + h)
    expect_identical(g$pvalue, b$pvalue)
    expect_equal(g$conf_int, b$conf_int, tolerance = 1e-7)
  }
})

test_that("an outcome that is not numeric or logical is refused by name", {
  dy <- make_dyadic(6, 6, seed = 62)
  expect_error(mwperm_dyadic(as.character(dy$y), dy$d, row = dy$i, col = dy$j),
               "`y` must be numeric", fixed = TRUE)
  expect_error(mwperm_dyadic(Sys.Date() + seq_len(nrow(dy)), dy$d, row = dy$i,
                             col = dy$j), "`y` must be numeric", fixed = TRUE)
})
