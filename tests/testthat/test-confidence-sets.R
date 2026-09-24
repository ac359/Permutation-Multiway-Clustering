## Confidence sets by test inversion: Procedure 1, step 3 of Guo, Toulis &
## Wang (2026) defines CI = {b : pval(b) > alpha}.
##
## The package computes that set in closed form for d = 1 (the p-value is a
## step function of b whose jumps are the roots of |v_k - W_k b| =
## |u_j - M_j b|; .ci_breakpoints()) and reports each component CLOSED (its
## end point is the bounding root, which may itself be rejected -- see D10).
## Duality is therefore checked just inside and just outside every reported
## end point, by refitting with that value as `beta_null` and the same seed:
## the refit runs exactly the groups the interval was inverted from.

fx <- make_dyadic(20, 20, beta = 0.3, seed = 61,
                  err = function(k) stats::rt(k, 3))
X <- ref_X(fx$x, nrow(fx))
R <- 5
seed <- 71
fit <- mwperm_dyadic(fx$y, fx$d, x = fx$x, row = fx$i, col = fx$j,
                     n_reps = R, seed = seed)
p_at <- function(b)
  mwperm_dyadic(fx$y, fx$d, x = fx$x, row = fx$i, col = fx$j, n_reps = R,
                seed = seed, beta_null = b, conf_int = FALSE)$pvalue

## The breakpoints of the aggregated step function, from the same groups the
## fit drew (rebuilt from the seed), so the probes can be placed strictly
## inside the open cells on either side of an end point.
roots <- local({
  preps <- lapply(seq_len(R), function(r) mwperm:::.ipt_prepare(
    fx$y, as.matrix(fx$d), X,
    ref_groups_crossed(list(fx$i, fx$j), fit$K, seed, r)))
  mwperm:::.ci_breakpoints(preps)
})
probe <- function(e) {
  others <- roots[abs(roots - e) > 1e-12 * (1 + abs(e))]
  min(abs(others - e)) / 10       # a tenth of the way to the next jump
}

test_that("the fit inverts the exact route", {
  expect_identical(fit$ci_method, "exact")
  expect_true(is.matrix(fit$conf_set))
  expect_identical(fit$conf_int, c(min(fit$conf_set), max(fit$conf_set)))
  ## every reported end point is one of the step function's jumps
  for (e in as.vector(fit$conf_set))
    expect_lt(min(abs(roots - e)), 1e-10 * (1 + abs(e)))
})

test_that("duality: pval > alpha just inside, <= alpha just outside", {
  a <- fit$alpha
  for (i in seq_len(nrow(fit$conf_set))) {
    lo <- fit$conf_set[i, 1]
    hi <- fit$conf_set[i, 2]
    expect_gt(p_at(lo + probe(lo)), a)
    expect_lte(p_at(lo - probe(lo)), a)
    expect_gt(p_at(hi - probe(hi)), a)
    expect_lte(p_at(hi + probe(hi)), a)
    expect_gt(p_at((lo + hi) / 2), a)                 # the interior
  }
})

test_that("the reported p-value is the one the set inverts", {
  ## 0 lies outside the set exactly when the test rejects beta = 0
  inside0 <- any(fit$conf_set[, 1] < 0 & 0 < fit$conf_set[, 2])
  expect_identical(inside0, fit$pvalue > fit$alpha)
  expect_identical(fit$pvalue, p_at(0))
})

test_that("reported end points are in {b : pval(b) > alpha} (paper's CI)", {
  skip_discrepancy("D10", paste("the exact route reports each component's",
             "CLOSURE, so an end point can be a rejected jump (e.g. the",
             "trade_dyadic anchor's lower end, p = 0.05 = alpha); the paper",
             "defines CI = {b : pval(b) > alpha}"))
  data(trade_dyadic, package = "mwperm", envir = environment())
  f <- with(trade_dyadic, mwperm_dyadic(log_trade, log_dist,
                                        x = cbind(log_gdp_i, log_gdp_j),
                                        row = importer, col = exporter,
                                        seed = 1))
  for (e in f$conf_int) {
    p <- with(trade_dyadic, mwperm_dyadic(
      log_trade, log_dist, x = cbind(log_gdp_i, log_gdp_j), row = importer,
      col = exporter, seed = 1, beta_null = e, conf_int = FALSE))$pvalue
    expect_gt(p, f$alpha)
  }
})

test_that("resolution guard: K + 1 < 1/alpha gives no set, with a note", {
  sm <- make_dyadic(15, 15, beta = 0.5, seed = 62)
  f <- mwperm_dyadic(sm$y, sm$d, x = sm$x, row = sm$i, col = sm$j,
                     n_reps = 3, seed = 1)
  expect_identical(f$K, 14L)
  expect_null(f$conf_int)
  expect_null(f$conf_set)
  expect_true(any(grepl("No 95% confidence interval", f$note, fixed = TRUE)))
  s <- NULL
  utils::capture.output(s <- summary(f))
  expect_true(is.na(s$ipt_ci_low) && is.na(s$ipt_ci_high))
  expect_error(confint(f), "at least 20 levels")
  ## the same design supports a 90% set, since 1/15 < 0.10
  f10 <- mwperm_dyadic(sm$y, sm$d, x = sm$x, row = sm$i, col = sm$j,
                       n_reps = 3, seed = 1, alpha = 0.10)
  expect_length(f10$conf_int, 2L)
  ## and median2 doubles the floor: K + 1 = 30 is too coarse for 95%
  m2 <- make_dyadic(30, 30, beta = 0.2, seed = 63)
  f2 <- mwperm_dyadic(m2$y, m2$d, x = m2$x, row = m2$i, col = m2$j,
                      n_reps = 3, seed = 1, aggregate = "median2")
  expect_null(f2$conf_int)
  expect_true(any(grepl("2/(K+1)", f2$note, fixed = TRUE)))
})

test_that("bisection and a user grid agree with the exact set", {
  ## Bisection (forced by a zero budget) stops at step * 1e-3; a grid is
  ## accurate to its spacing and reports attained, accepted grid points.
  old <- options(mwperm.ci_exact_budget = 0)
  on.exit(options(old))
  fb <- mwperm_dyadic(fx$y, fx$d, x = fx$x, row = fx$i, col = fx$j,
                      n_reps = R, seed = seed)
  options(old)
  expect_identical(fb$ci_method, "bisection")
  tol <- fit$se_naive * 1e-3
  expect_lt(max(abs(fb$conf_int - fit$conf_int)), tol)

  h <- 0.002
  grid <- seq(fit$conf_int[1] - 0.2, fit$conf_int[2] + 0.2, by = h)
  fg <- mwperm_dyadic(fx$y, fx$d, x = fx$x, row = fx$i, col = fx$j,
                      n_reps = R, seed = seed, grid = grid)
  expect_identical(fg$ci_method, "grid")
  expect_true(all(abs(fg$conf_int - fit$conf_int) <= h + 1e-12))
  expect_gte(fg$conf_int[1], fit$conf_int[1])   # grid points are attained,
  expect_lte(fg$conf_int[2], fit$conf_int[2])   # so they lie inside the hull
  expect_true(any(grepl("accurate to the grid spacing", fg$note)))
})

test_that("d > 1: retained region points are accepted by the joint test", {
  D2 <- cbind(d1 = fx$d, d2 = with_seed(64, stats::rnorm(nrow(fx))))
  y <- fx$y + 0.2 * D2[, 2]
  g <- list(seq(0, 0.6, by = 0.1), seq(-0.1, 0.5, by = 0.1))
  f <- mwperm_dyadic(y, D2, x = fx$x, row = fx$i, col = fx$j, n_reps = 3,
                     seed = 5, grid = g)
  expect_null(f$conf_int)
  pts <- f$conf_region
  expect_gt(nrow(pts), 0L)
  pj <- function(b) mwperm_dyadic(y, D2, x = fx$x, row = fx$i, col = fx$j,
                                  n_reps = 3, seed = 5, beta_null = b,
                                  conf_int = FALSE)$pvalue
  for (i in unique(round(seq(1, nrow(pts), length.out = 4))))
    expect_gt(pj(pts[i, ]), f$alpha)
  ## a grid point outside the retained set is rejected
  all_pts <- as.matrix(expand.grid(g))
  out <- all_pts[!paste(all_pts[, 1], all_pts[, 2]) %in%
                   paste(pts[, 1], pts[, 2]), , drop = FALSE]
  expect_lte(pj(out[1, ]), f$alpha)
  expect_equal(unname(f$conf_box[1, ]), unname(apply(pts, 2, min)))
})
