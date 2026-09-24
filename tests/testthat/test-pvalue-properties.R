## Properties of the randomization p-value of Eq. (10) (Guo, Toulis & Wang
## 2026) as the front ends report it.
##   * each repetition's p-value is (1 + #{k : min_j a_j <= b_k}) / (K + 1),
##     so it lies on the grid {1, ..., K + 1} / (K + 1);
##   * the reported p-value is the median over repetitions (Remark 1) -- with
##     an even n_reps that is the mean of the two central values, which can
##     fall between grid points (the draft says so in Section 2.5);
##   * ties count toward the p-value: the indicator is <=, not <;
##   * if D lies in col(X), then V_k' D = 0 for every k, so a_k = b_k = 0 and
##     Eq. (10) gives p = 1 exactly.
## Complements tests/lower-level-tests/test-pvalue.R, which checks the same
## grid and tie rules on hand-built internal objects rather than on fits.

on_grid <- function(p, order) {
  m <- p * order
  all(abs(m - round(m)) < 1e-9) && all(p >= 1 / order - 1e-12) &&
    all(p <= 1 + 1e-12)
}

test_that("every per-repetition p-value is a multiple of 1/(K+1) in [1/(K+1), 1]", {
  dy <- make_dyadic(12, 9, beta = 0.2, seed = 51)
  pn <- make_panel(9, 8, T = 3, beta = 0.2, seed = 52)
  tw <- make_threeway(6, 7, 8, beta = 0.2, seed = 53)
  ly <- make_layout(4, 4, sizes = 7:9, beta = 0.3, seed = 54)
  fits <- list(
    mwperm_dyadic(dy$y, dy$d, x = dy$x, row = dy$i, col = dy$j,
                  n_reps = 7, seed = 1, conf_int = FALSE),
    mwperm_panel(pn$y, pn$d, x = pn$x, row = pn$i, col = pn$j, time = pn$t,
                 n_reps = 7, seed = 1, conf_int = FALSE),
    mwperm_threeway(tw$y, tw$d, x = tw$x, id1 = tw$i, id2 = tw$j, id3 = tw$l,
                    n_reps = 7, seed = 1, conf_int = FALSE),
    mwperm_layout(ly$y, ly$d, x = ly$x, row = ly$i, col = ly$j, rep = ly$l,
                  n_reps = 7, seed = 1, conf_int = FALSE),
    with(subset(dy, i != j),
         mwperm_missing(y, d, x = x, row = i, col = j, min_block = 3,
                        n_reps = 7, seed = 1, conf_int = FALSE)))
  for (f in fits) {
    expect_true(on_grid(f$pvalues_rep, f$K + 1))
    expect_identical(f$n_perm, f$K + 1L)
    expect_equal(f$resolution, 1 / (f$K + 1))
  }
  ## the sign-flip group: order 2^(n_flip - 1)
  h <- mwperm_dyadic_het(dy$y, dy$d, x = dy$x, row = dy$i, col = dy$j,
                         n_flip = 5, n_reps = 7, seed = 1, conf_int = FALSE)
  expect_identical(h$n_perm, 16L)
  expect_true(on_grid(h$pvalues_rep, 16))
})

test_that("the reported p-value is the median over repetitions, odd or even", {
  dy <- make_dyadic(10, 10, beta = 0.12, seed = 55)
  off_grid_seen <- FALSE
  for (R in c(1, 2, 4, 5, 10)) for (s in 1:6) {
    f <- mwperm_dyadic(dy$y, dy$d, x = dy$x, row = dy$i, col = dy$j,
                       n_reps = R, seed = s, conf_int = FALSE)
    expect_identical(f$pvalue, stats::median(f$pvalues_rep))
    if (!on_grid(f$pvalue, f$K + 1)) {
      off_grid_seen <- TRUE
      expect_identical(R %% 2, 0)       # only an even n_reps can leave the grid
      expect_true(any(grepl("can fall between grid points",
                            utils::capture.output(print(f)), fixed = TRUE)))
    }
  }
  expect_true(off_grid_seen)            # the even-median case does occur
})

test_that("aggregate = 'median2' reports min(1, 2 x median)", {
  dy <- make_dyadic(10, 10, beta = 0.12, seed = 55)
  f <- mwperm_dyadic(dy$y, dy$d, x = dy$x, row = dy$i, col = dy$j,
                     n_reps = 5, seed = 3, conf_int = FALSE,
                     aggregate = "median2")
  expect_identical(f$pvalue, min(1, 2 * stats::median(f$pvalues_rep)))
  expect_equal(f$p_floor, 2 / (f$K + 1))
})

test_that("ties count toward the p-value (<=, not <)", {
  ## y = 0 exactly: every residual inner product is an exact zero, so
  ## a_k = b_k = 0 for all k and every indicator in Eq. (10) is a tie.
  ## With <= the p-value is (1 + K)/(K + 1) = 1; with < it would be 1/(K+1),
  ## a rejection of H0 on data with no signal at all.
  dy <- make_dyadic(10, 10, seed = 56)
  f <- mwperm_dyadic(numeric(nrow(dy)), dy$d, x = dy$x, row = dy$i,
                     col = dy$j, n_reps = 3, seed = 2, conf_int = FALSE)
  expect_identical(f$pvalues_rep, c(1, 1, 1))
  ## A layout whose d is constant within every cell: every within-cell
  ## permutation leaves the residualized d unchanged, so b_k = a_k exactly
  ## for every k, the minimum a_j is matched by its own b_j, and p = 1.
  ly <- make_layout(4, 4, sizes = 6:8, seed = 57, d_cell_constant = TRUE)
  expect_warning(
    fl <- mwperm_layout(ly$y, ly$d, x = ly$x, row = ly$i, col = ly$j,
                        rep = ly$l, n_reps = 3, seed = 2, conf_int = FALSE),
    "constant within every cell")
  expect_identical(fl$pvalues_rep, c(1, 1, 1))
})

test_that("D in col(X): p = 1 exactly, as Eq. (10) gives (and a warning)", {
  ## GTW: V_k' X = 0 and D = X c imply V_k' D = 0, hence a_k = b_k = 0 and
  ## pval = 1; the confidence set is then the whole line, {b : 1 > alpha}.
  ## The package matches the paper and warns that beta is unidentified.
  dy <- make_dyadic(20, 20, seed = 58)
  d_in_x <- 3 - 2 * dy$x
  expect_warning(
    f <- mwperm_dyadic(dy$y, d_in_x, x = dy$x, row = dy$i, col = dy$j,
                       n_reps = 3, seed = 2),
    "no variation after partialling out")
  expect_identical(f$pvalues_rep, c(1, 1, 1))
  expect_identical(f$conf_int, c(-Inf, Inf))
  ## a constant d is the same case: it lies in the span of the intercept
  expect_warning(
    f2 <- mwperm_dyadic(dy$y, rep(5, nrow(dy)), x = dy$x, row = dy$i,
                        col = dy$j, n_reps = 2, seed = 2, conf_int = FALSE),
    "no variation")
  expect_identical(f2$pvalues_rep, c(1, 1))
})
