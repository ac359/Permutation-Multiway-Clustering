## The "mwperm" object and its S3 methods as the JSS draft documents them
## (Section 3.5): the principal slots exist with the right types; summary()
## is one row per coefficient; confint() carries the "2.5 %"/"97.5 %" columns
## and its "method" attribute; print() is pinned by a snapshot; plot() draws
## for one and for many repetitions. tests/test-methods.R (base R) pins the
## label contract in more detail; this file checks the draft's statements.

fx <- make_dyadic(20, 20, beta = 0.3, seed = 111)
fit <- mwperm_dyadic(fx$y, fx$d, x = fx$x, row = fx$i, col = fx$j,
                     n_reps = 5, seed = 1)

test_that("the slots named in draft Sec. 3.5 exist with the right types", {
  expect_s3_class(fit, "mwperm")
  expect_type(fit$estimate, "double")
  expect_named(fit$estimate, "fx$d")
  expect_true(is.numeric(fit$pvalue) && length(fit$pvalue) == 1L)
  expect_true(fit$pvalue > 0 && fit$pvalue <= 1)
  expect_type(fit$pvalues_rep, "double")
  expect_length(fit$pvalues_rep, 5L)
  expect_identical(fit$aggregate, "median")
  expect_true(is.matrix(fit$conf_set) && ncol(fit$conf_set) == 2L)
  expect_type(fit$conf_int, "double")
  expect_length(fit$conf_int, 2L)
  expect_identical(fit$alpha, 0.05)
  expect_identical(fit$conf_level, 0.95)
  expect_type(fit$note, "character")
  ## the fields the methods read
  for (nm in c("se_naive", "ci_method", "K", "n_perm", "n_reps", "type",
               "d_names", "n_obs", "n_clusters", "resolution", "p_floor",
               "beta_null", "call"))
    expect_true(nm %in% names(fit), info = nm)
  ## biclique designs add cells_used; the sign-flip test adds n_flip
  md <- subset(fx, i != j)
  fm <- mwperm_missing(md$y, md$d, x = md$x, row = md$i, col = md$j,
                       n_reps = 1, seed = 1, conf_int = FALSE)
  expect_type(fm$cells_used, "integer")
  expect_lt(nobs(fm), nrow(md))           # the test used fewer cells
  fh <- mwperm_dyadic_het(fx$y, fx$d, row = fx$i, col = fx$j, n_flip = 4,
                          n_reps = 1, seed = 1, conf_int = FALSE)
  expect_identical(fh$n_flip, 4L)
})

test_that("summary() returns one row per coefficient", {
  s <- NULL
  utils::capture.output(s <- summary(fit))
  expect_s3_class(s, "data.frame")
  expect_identical(nrow(s), 1L)
  expect_named(s, c("term", "ols_estimate", "ols_se_naive", "ipt_ci_low",
                    "ipt_ci_high", "p_value"))
  expect_identical(s$ipt_ci_low, fit$conf_int[1])
  D2 <- cbind(a = fx$d, b = with_seed(112, stats::rnorm(nrow(fx))))
  f2 <- mwperm_dyadic(fx$y, D2, x = fx$x, row = fx$i, col = fx$j,
                      n_reps = 1, seed = 1, conf_int = FALSE)
  utils::capture.output(s2 <- summary(f2))
  expect_identical(nrow(s2), 2L)
  expect_identical(s2$term, c("a", "b"))
})

test_that("confint() has percentile columns and names its method", {
  ci <- confint(fit)
  expect_identical(colnames(ci), c("2.5 %", "97.5 %"))
  expect_identical(attr(ci, "method"), "IPT (inverted permutation test)")
  expect_identical(unname(ci[1, ]), fit$conf_int)
  f90 <- mwperm_dyadic(fx$y, fx$d, x = fx$x, row = fx$i, col = fx$j,
                       n_reps = 1, seed = 1, alpha = 0.1)
  expect_identical(colnames(confint(f90)), c("5 %", "95 %"))
  expect_error(confint(fit, level = 0.9), "refitting")
})

test_that("coef() and nobs() report the OLS estimate and the cells used", {
  expect_identical(coef(fit), fit$estimate)
  expect_identical(nobs(fit), nrow(fx))
  ols <- stats::lm(y ~ d + x, data = fx)
  expect_equal(unname(coef(fit)), unname(stats::coef(ols)["d"]),
               tolerance = 1e-10)
})

test_that("print() output is stable (snapshot)", {
  local_reproducible_output(width = 80)
  expect_snapshot(print(fit))
  md <- subset(fx, i != j)
  expect_snapshot(print(mwperm_missing(md$y, md$d, x = md$x, row = md$i,
                                       col = md$j, n_reps = 3, seed = 1)))
  expect_snapshot(print(mwperm_dyadic_het(fx$y, fx$d, x = fx$x, row = fx$i,
                                          col = fx$j, n_flip = 6,
                                          n_reps = 3, seed = 1)))
})

test_that("plot() draws for one and for many repetitions", {
  f1 <- mwperm_dyadic(fx$y, fx$d, x = fx$x, row = fx$i, col = fx$j,
                      n_reps = 1, seed = 1)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  settable <- c("mar", "mfrow", "mfcol", "oma", "mgp", "las", "cex", "tcl",
                "xpd", "bty", "family", "font.main")
  op <- graphics::par(settable)
  expect_no_error(plot(f1))
  expect_no_error(plot(fit))
  expect_no_error(plot(fit, type = "stability"))
  expect_no_error(plot(fit, type = "all"))
  expect_identical(graphics::par(settable), op)   # par() restored
})
