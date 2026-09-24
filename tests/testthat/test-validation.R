## Input validation, and the equivalence of the calling interfaces.
##
## Theorem 1 of Guo, Toulis & Wang (2026) is stated for p < N/2 (Procedure 1
## needs an (N - 2p)-column V_k); the engine refuses N <= 2p. Every other
## check here guards the data contract: bad input must fail early, naming the
## argument, rather than inside a factorization. tests/test-validation.R
## (base R) covers the scalar-argument checks; this file covers the data.

dy <- make_dyadic(8, 8, seed = 121)

test_that("mismatched lengths name the offending argument", {
  expect_error(mwperm_dyadic(dy$y, dy$d, row = dy$i[-1], col = dy$j),
               "`row` must have the same length as `y`")
  expect_error(mwperm_dyadic(dy$y, dy$d[-1], row = dy$i, col = dy$j),
               "`d` must have the same length|same number of rows|length")
  expect_error(mwperm_dyadic(dy$y, dy$d, x = dy$x[-1], row = dy$i,
                             col = dy$j),
               "`x` must have the same number of rows as `y`")
})

test_that("NA / non-finite values are refused, by argument", {
  y <- dy$y
  y[3] <- NA
  expect_error(mwperm_dyadic(y, dy$d, row = dy$i, col = dy$j, seed = 1),
               "`y` contains missing or non-finite")
  d <- dy$d
  d[5] <- Inf
  expect_error(mwperm_dyadic(dy$y, d, row = dy$i, col = dy$j, seed = 1),
               "`d` contains missing or non-finite")
  x <- dy$x
  x[2] <- NaN
  expect_error(mwperm_dyadic(dy$y, dy$d, x = x, row = dy$i, col = dy$j,
                             seed = 1), "`x` contains missing or non-finite")
  i <- dy$i
  i[1] <- NA
  expect_error(mwperm_dyadic(dy$y, dy$d, row = i, col = dy$j, seed = 1),
               "`row` contains missing values")
})

test_that("a non-numeric d or y is refused", {
  expect_error(mwperm_dyadic(dy$y, as.character(dy$d), row = dy$i,
                             col = dy$j, seed = 1), "`d` must be numeric")
  expect_error(mwperm_dyadic(factor(round(dy$y)), dy$d, row = dy$i,
                             col = dy$j, seed = 1), "`y` is a factor")
})

test_that("fewer than two clusters on a dimension is refused", {
  one <- subset(dy, i == 1)
  expect_error(mwperm_dyadic(one$y, one$d, row = one$i, col = one$j,
                             seed = 1), "need >= 2 in each dimension")
  expect_error(mwperm_check(index = one[c("i", "j")]),
               "Fewer than 2 effective clustering dimensions")
})

test_that("duplicate (i, j) cells without `rep` are refused, or read as a layout", {
  dup <- rbind(dy, dy[1:5, ])
  expect_error(mwperm_dyadic(dup$y, dup$d, row = dup$i, col = dup$j,
                             seed = 1), "one observation per \\(row, col\\)")
  expect_error(mwperm_missing(dup$y, dup$d, row = dup$i, col = dup$j,
                              seed = 1), "at most one observation")
  ## the dispatcher reads repeated cells as within-cell replication, and says
  ## so with a warning, because that is an assumption about the data
  expect_warning(
    expect_identical(suppressMessages(mwperm_check(
      index = dup[c("i", "j")]))$design, "layout"),
    NA)
  lay <- make_layout(4, 4, sizes = 6:7, seed = 122)
  expect_warning(suppressMessages(mwperm(lay$y, lay$d,
                                         index = lay[c("i", "j")],
                                         n_reps = 1, seed = 1,
                                         conf_int = FALSE)),
                 "treated as exchangeable within-cell replication")
})

test_that("p >= N/2 is refused (Theorem 1 needs p < N/2)", {
  sm <- make_dyadic(4, 4, seed = 123)                 # N = 16
  x7 <- with_seed(124, matrix(stats::rnorm(16 * 7), 16, 7))
  expect_error(mwperm_dyadic(sm$y, sm$d, x = x7, row = sm$i, col = sm$j,
                             seed = 1), "Need N > 2p")   # p = 8, N = 2p
})

test_that("the vector and the column-name interfaces agree", {
  fx <- make_dyadic(20, 20, beta = 0.2, seed = 125)
  v <- suppressMessages(mwperm(fx$y, fx$d, x = fx$x, index = fx[c("i", "j")],
                               n_reps = 3, seed = 1))
  n <- suppressMessages(mwperm(y = "y", d = "d", x = "x",
                               index = c("i", "j"), data = fx, n_reps = 3,
                               seed = 1))
  f <- suppressMessages(mwperm_formula(y ~ d | x, data = fx,
                                       index = c("i", "j"), n_reps = 3,
                                       seed = 1))
  for (o in list(n, f)) {
    expect_identical(o$pvalues_rep, v$pvalues_rep)
    expect_identical(o$conf_int, v$conf_int)
    expect_identical(unname(o$estimate), unname(v$estimate))
  }
})

test_that("factor, character and integer cluster ids agree when they sort alike", {
  ## The dense coding follows factor() level order, so ids of any type give
  ## the same groups -- and the same seeded result -- when their levels sort
  ## in the same order. Character ids sort by collation ("10" < "2"), which
  ## is a different, equally valid relabelling: the draft (Sec. 3.8) and
  ## ?build_perm_set say so.
  fx <- make_dyadic(12, 12, beta = 0.2, seed = 126)
  base <- mwperm_dyadic(fx$y, fx$d, row = fx$i, col = fx$j, n_reps = 3,
                        seed = 1, conf_int = FALSE)
  variants <- list(
    factor = list(factor(fx$i), factor(fx$j)),
    padded = list(sprintf("r%02d", fx$i), sprintf("c%02d", fx$j)),
    double = list(as.double(fx$i), as.double(fx$j)),
    labels = list(factor(fx$i, labels = LETTERS[1:12]),
                  factor(fx$j, labels = letters[1:12])))
  for (v in variants) {
    f <- mwperm_dyadic(fx$y, fx$d, row = v[[1]], col = v[[2]], n_reps = 3,
                       seed = 1, conf_int = FALSE)
    expect_identical(f$pvalues_rep, base$pvalues_rep)
  }
  ## unpadded character ids: a different relabelling, same estimate
  ch <- mwperm_dyadic(fx$y, fx$d, row = as.character(fx$i),
                      col = as.character(fx$j), n_reps = 3, seed = 1,
                      conf_int = FALSE)
  expect_identical(unname(ch$estimate), unname(base$estimate))
  expect_identical(ch$K, base$K)
})
