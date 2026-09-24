## Procedure 1 and Eq. (10) of Guo, Toulis & Wang (2026), checked against the
## naive implementation in helper-reference.R.
##
## Two levels, for the same permutation set on both sides:
##   * statistics: the package's a_k and b_k (read from its cached cross
##     products, .ipt_prepare() / .ipt_eval()) against the naive
##     ||D' V_k V_k' y|| and ||D' V_k V_k' y_k||, relative tolerance 1e-8;
##   * p-values: every entry of a front end's `pvalues_rep` against the naive
##     Eq. (10) on the group that repetition drew, rebuilt from the seed.
## Near-ties: an indicator 1{min_j a_j <= b_k} decided within rounding could
## flip between two correct implementations, so every p-value comparison
## first asserts that the reference's nearest indicator is more than 1e-6
## from a tie (the designs below were chosen with that margin), and the
## statistics themselves are compared directly.

## The package's statistics for a given list of gather vectors.
pkg_stats <- function(y, D, X, perms, b = 0) {
  prep <- mwperm:::.ipt_prepare(y, as.matrix(D), X, perms, need_perm_D = TRUE)
  mwperm:::.ipt_eval(prep, b)
}

expect_matches_reference <- function(fit, y, D, X, groups_of, b = 0) {
  n_reps <- length(fit$pvalues_rep)
  ref <- ref_reps(y, D, X, groups_of, n_reps, b = b)
  for (r in seq_len(n_reps)) {
    perms <- groups_of(r)
    pk <- pkg_stats(y, D, X, perms, b = b)
    expect_equal(pk$a, ref[[r]]$a, tolerance = 1e-8)
    expect_equal(pk$b, ref[[r]]$b, tolerance = 1e-8)
    expect_gt(ref[[r]]$gap, 1e-6)
    expect_identical(fit$pvalues_rep[r], ref[[r]]$pvalue)
  }
  invisible(ref)
}

test_that("dyadic: a_k, b_k and every per-repetition p-value match Eq. (10)", {
  dat <- make_dyadic(12, 10, beta = 0.15, seed = 11)
  K <- 7L
  fit <- mwperm_dyadic(dat$y, dat$d, x = dat$x, row = dat$i, col = dat$j,
                       K = K, n_reps = 4, seed = 21, conf_int = FALSE)
  X <- ref_X(dat$x, nrow(dat))
  groups <- function(r) ref_groups_crossed(list(dat$i, dat$j), K, 21, r)
  expect_matches_reference(fit, dat$y, dat$d, X, groups)
})

test_that("the package builds y_{pi_k, sigma_k} as GTW define it", {
  ## GTW Section 3.1: the entry of y_{pi,sigma} for cell (i, j) is
  ## y_{pi(i) sigma(j)}. The package's gather vectors (.build_obs_perms, the
  ## mixed-radix cell codes) must equal the reference's key-matched ones.
  dat <- make_dyadic(9, 7, seed = 3)
  ri <- dense_id(dat$i)
  ci <- dense_id(dat$j)
  Gr <- build_perm_set(9, 6, seed = 5)
  Gc <- build_perm_set(7, 6, seed = 6)
  pkg <- mwperm:::.build_obs_perms(cbind(ri, ci), list(Gr, Gc))
  ref <- lapply(1:7, function(k) ref_gather(list(ri, ci),
                                            list(Gr[[k]], Gc[[k]])))
  expect_identical(pkg, ref)
})

test_that("the lexicographic stacking of Eq. (7) is immaterial", {
  ## GTW stack cell (i, j) in row (i - 1) n + j. The package keys cells by
  ## code instead, so a fit on shuffled rows must give the same p-values as
  ## the reference run on the lexicographic stacking.
  dat <- make_dyadic(10, 10, beta = 0.2, seed = 12)
  lex <- dat[order(dat$i, dat$j), ]
  shuf <- lex[with_seed(99, sample.int(nrow(lex))), ]
  fit <- mwperm_dyadic(shuf$y, shuf$d, x = shuf$x, row = shuf$i,
                       col = shuf$j, K = 9, n_reps = 3, seed = 4,
                       conf_int = FALSE)
  groups <- function(r) ref_groups_crossed(list(lex$i, lex$j), 9, 4, r)
  ref <- ref_reps(lex$y, lex$d, ref_X(lex$x, nrow(lex)), groups, 3)
  expect_identical(fit$pvalues_rep, vapply(ref, `[[`, 0, "pvalue"))
})

test_that("a non-zero null is Procedure 1 on y - D b", {
  dat <- make_dyadic(11, 11, beta = 0.4, seed = 13)
  fit <- mwperm_dyadic(dat$y, dat$d, x = dat$x, row = dat$i, col = dat$j,
                       K = 10, n_reps = 3, seed = 8, beta_null = 0.35,
                       conf_int = FALSE)
  groups <- function(r) ref_groups_crossed(list(dat$i, dat$j), 10, 8, r)
  expect_matches_reference(fit, dat$y, dat$d, ref_X(dat$x, nrow(dat)),
                           groups, b = 0.35)
})

test_that("d > 1: the joint statistic ||D' V_k V_k' y|| matches", {
  dat <- make_dyadic(10, 9, seed = 14)
  D2 <- cbind(d1 = dat$d, d2 = with_seed(15, stats::rnorm(nrow(dat))))
  y <- dat$y + 0.3 * D2[, 2]
  fit <- mwperm_dyadic(y, D2, x = dat$x, row = dat$i, col = dat$j, K = 8,
                       n_reps = 3, seed = 16, conf_int = FALSE)
  groups <- function(r) ref_groups_crossed(list(dat$i, dat$j), 8, 16, r)
  expect_matches_reference(fit, y, D2, ref_X(dat$x, nrow(dat)), groups)
  ## ... and at a joint null vector
  fit_b <- mwperm_dyadic(y, D2, x = dat$x, row = dat$i, col = dat$j, K = 8,
                         n_reps = 3, seed = 16, conf_int = FALSE,
                         beta_null = c(0.1, 0.25))
  expect_matches_reference(fit_b, y, D2, ref_X(dat$x, nrow(dat)), groups,
                           b = c(0.1, 0.25))
})

test_that("rank-deficient [X | X_k]: row and column fixed effects", {
  ## Row dummies are permuted into row dummies, so col(X_k) and col(X) share
  ## every fixed-effect direction and the stack is far from rank 2p. The
  ## reference uses the numerical rank of the stack; so must the package.
  dat <- make_dyadic(10, 10, beta = 0.3, seed = 17)
  fe <- cbind(stats::model.matrix(~ factor(dat$i))[, -1],
              stats::model.matrix(~ factor(dat$j))[, -1])
  x <- cbind(dat$x, fe)
  X <- ref_X(x, nrow(dat))
  M <- cbind(X, X[ref_groups_crossed(list(dat$i, dat$j), 9, 5, 1)[[2]], ])
  expect_lt(qr(M)$rank, ncol(M))          # genuinely rank deficient
  fit <- mwperm_dyadic(dat$y, dat$d, x = x, row = dat$i, col = dat$j,
                       K = 9, n_reps = 3, seed = 5, conf_int = FALSE)
  groups <- function(r) ref_groups_crossed(list(dat$i, dat$j), 9, 5, r)
  expect_matches_reference(fit, dat$y, dat$d, X, groups)
})

test_that("the minimum is over the K non-identity elements only (Eq. 10)", {
  ## Eq. (10) takes min_{1 <= j <= K} a_j; including the identity (whose
  ## stacked design [X | X] has rank p, not 2p) would lower the minimum.
  ## The package's a-vector has exactly K entries.
  dat <- make_dyadic(8, 8, seed = 18)
  perms <- ref_groups_crossed(list(dat$i, dat$j), 7, 3, 1)
  pk <- pkg_stats(dat$y, dat$d, ref_X(dat$x, nrow(dat)), perms)
  expect_length(pk$a, 7L)
  expect_length(pk$b, 7L)
})
