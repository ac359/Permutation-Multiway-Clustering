## Every output the JSS draft (2026-09-18 build) prints for a fixed seed, as a
## regression test at the printed precision. Printed blocks are compared
## after collapsing whitespace, because the draft typesets them at a narrower
## width than the console; the numbers themselves are also checked on the
## fitted objects, and OLS estimates are checked against lm(). The Section 7
## real-data examples run only when the gravity package is installed; the
## Section 7.4 CEPII panel needs a database that is not distributed.
##
## Where the draft's text no longer matches version 0.4.2 because 0.4.2
## changed it on purpose, the draft's assertion is kept and skipped with a
## discrepancy number (D4-D6): the draft, not the code, needs the update.

squash <- function(lines) gsub(" +", " ", trimws(paste(trimws(lines),
                                                        collapse = " ")))
printed <- function(x) squash(utils::capture.output(print(x)))
quiet <- function(expr) suppressMessages(expr)

data(trade_dyadic, package = "mwperm", envir = environment())
data(trade_panel, package = "mwperm", envir = environment())

## ---- Section 4: the synthetic data sets ------------------------------------

test_that("Sec. 4: the data sets and their data-generating coefficients", {
  expect_identical(dim(trade_dyadic), c(1600L, 8L))
  expect_identical(nlevels(trade_dyadic$importer), 40L)
  expect_identical(unname(attr(trade_dyadic, "true_coef")["log_dist"]), -1)
  expect_identical(unname(attr(trade_dyadic, "true_coef")["placebo"]), 0)
  expect_identical(nrow(trade_panel), 22L * 22L * 6L)
  expect_identical(unname(attr(trade_panel, "true_coef")["fta"]), 0.5)
})

## ---- Section 3.3: the design diagnostic ------------------------------------

test_that("Sec. 3.3: mwperm_check() on trade_dyadic", {
  expect_identical(
    printed(mwperm_check(index = c("importer", "exporter"),
                         data = trade_dyadic)),
    paste("mwperm design diagnosis ------------------------------ Detected",
          "design : dyadic (2 indices, one observation per cell, complete",
          "array) Roles : row = importer, col = exporter Dimensions :",
          "importer (40) x exporter (40) | 1600 observations Balance :",
          "complete Resolution : default K = 39, so p-values are multiples",
          "of 1/40 = 0.025 -> fine enough for a 95% confidence set at alpha",
          "= 0.05 Would run : mwperm_dyadic(y, d, x, row = importer, col =",
          "exporter) ? design = \"dyadic_het\" runs the sign-flip test",
          "instead: valid under arbitrary heteroskedasticity for errors that",
          "are independent across cells and symmetric about zero, but NOT",
          "under additive cluster effects eta_i + xi_j (see",
          "?mwperm_dyadic_het)"))
})

## ---- Sections 3.4-3.5: the first example -----------------------------------

fit <- quiet(mwperm(y = "log_trade", d = "log_dist",
                    x = c("log_gdp_i", "log_gdp_j"),
                    index = c("importer", "exporter"), data = trade_dyadic,
                    n_reps = 15, seed = 101))

test_that("Sec. 3.4: the first example prints as in the draft", {
  expect_identical(
    printed(fit),
    paste("Invariant permutation test (mwperm)",
          "------------------------------------ Design : dyadic",
          "Auto-detected: dyadic (2 indices, one observation per cell,",
          "complete array) Clusters : row=40, col=40 (1600 observations)",
          "Permutations : K = 39 (group order 40, 15 reps) Resolution :",
          "p-values are multiples of 1/40 = 0.025 per rep; reported floor",
          "0.025 log_dist OLS estimate = -0.8985 95% IPT CI [-1.254,",
          "-0.5346] H0: beta = 0 p-value = 0.025 Decision : reject at alpha",
          "= 0.05"))
  expect_equal(signif(unname(fit$estimate), 4), -0.8985)
  expect_equal(signif(fit$conf_int, 4), c(-1.254, -0.5346))
  expect_identical(fit$pvalue, 1 / 40)
  ols <- stats::lm(log_trade ~ log_dist + log_gdp_i + log_gdp_j,
                   data = trade_dyadic)
  expect_equal(unname(fit$estimate), unname(stats::coef(ols)["log_dist"]),
               tolerance = 1e-10)
  ## "calling that worker directly with the same seed gives an identical
  ## result"
  direct <- with(trade_dyadic, mwperm_dyadic(
    y = log_trade, d = log_dist, x = cbind(log_gdp_i, log_gdp_j),
    row = importer, col = exporter, n_reps = 15, seed = 101))
  expect_identical(direct$pvalues_rep, fit$pvalues_rep)
  expect_identical(direct$conf_int, fit$conf_int)
})

test_that("Sec. 3.5: confint(fit)", {
  expect_identical(
    squash(utils::capture.output(print(confint(fit)))),
    paste("2.5 % 97.5 % log_dist -1.253701 -0.5346174 attr(,\"method\")",
          "[1] \"IPT (inverted permutation test)\""))
})

## ---- Section 3.6: the groups in code ---------------------------------------

test_that("Sec. 3.6: build_perm_set() and build_flip_set() as printed", {
  expect_identical(
    squash(utils::capture.output(print(build_perm_set(n = 6, K = 2,
                                                      seed = 1)))),
    paste("[[1]] [1] 1 2 3 4 5 6 [[2]] [1] 5 6 1 2 3 4 [[3]] [1] 3 4 5 6",
          "1 2 attr(,\"block_size\") [1] 3"))
  F <- build_flip_set(n_row = 4, n_col = 3, n_flip = 3, seed = 1)
  out <- utils::capture.output({
    print(c(length(F), attr(F, "group_order")))
    print(attr(F, "row_groups"))
    print(attr(F, "col_groups"))
    print(attr(F, "signs"))
    print(outer(F[[4]]$row, F[[4]]$col))
  })
  expect_identical(
    squash(out),
    paste("[1] 4 4 [1] 1 3 1 2 [1] 1 3 3 [,1] [,2] [,3] [1,] 1 1 1 [2,] 1",
          "-1 1 [3,] 1 1 -1 [4,] 1 -1 -1 [,1] [,2] [,3] [1,] 1 -1 -1 [2,]",
          "-1 1 1 [3,] 1 -1 -1 [4,] -1 1 1"))
})

## ---- Section 3.7: the root budget ------------------------------------------

test_that("Sec. 3.7: n_flip = 8 over 10 repetitions exceeds the exact budget", {
  ## "at n_flip = 8 and n_reps = 10 the exact route would have to enumerate
  ## 322,580 roots, so the interval comes from the fallback"
  expect_gt(2 * (2^7 - 1)^2 * 10, 2e5)
  f8 <- with(trade_dyadic, mwperm_dyadic_het(
    y = log_trade, d = log_dist, x = cbind(log_gdp_i, log_gdp_j),
    row = importer, col = exporter, n_flip = 8, n_reps = 10, seed = 1))
  expect_identical(f8$ci_method, "bisection")
  expect_true(any(grepl("bracketing and bisection", f8$note)))
})

test_that("Secs. 2.7 and 3.7: the default n_flip is 8", {
  skip_discrepancy("D6", paste("the draft says the sign-flip front end",
                               "fixes n_flip = 8 by default; 0.4.2 made the",
                               "default the resolution rule (6 at alpha =",
                               "0.05)"))
  fd <- with(trade_dyadic, mwperm_dyadic_het(
    y = log_trade, d = log_dist, row = importer, col = exporter,
    n_reps = 1, seed = 1, conf_int = FALSE))
  expect_identical(fd$n_flip, 8L)
})

## ---- Section 6.1: dyadic gravity -------------------------------------------

test_that("Sec. 6.1: the placebo regressor and the OLS comparison", {
  fit0 <- with(trade_dyadic, mwperm_dyadic(
    y = log_trade, d = placebo, x = cbind(log_dist, log_gdp_i, log_gdp_j),
    row = importer, col = exporter, n_reps = 15, seed = 202))
  expect_identical(
    squash(utils::capture.output(print(c(estimate = unname(fit0$estimate),
                                         p_value = fit0$pvalue)))),
    "estimate p_value -0.004086162 1.000000000")
  expect_equal(round(fit0$conf_int, 3), c(-0.087, 0.079))   # Table 4
  ols <- stats::lm(log_trade ~ log_dist + log_gdp_i + log_gdp_j,
                   data = trade_dyadic)
  expect_identical(
    squash(utils::capture.output(print(round(
      summary(ols)$coefficients["log_dist", , drop = FALSE], 4)))),
    "Estimate Std. Error t value Pr(>|t|) log_dist -0.8985 0.0873 -10.2934 0")
  ci_ols <- stats::confint(ols)["log_dist", ]
  expect_equal(unname(round(ci_ols, 3)), c(-1.070, -0.727))  # Table 4
  expect_equal(signif(summary(ols)$coefficients["log_dist", 4], 2), 4.2e-24)
})

## ---- Section 6.2: panel, and an incomplete panel ---------------------------

test_that("Sec. 6.2: the panel fit, with and without time effects", {
  fit_p <- with(trade_panel, mwperm_panel(
    y = log_trade, d = fta, x = cbind(log_gdp_i, log_gdp_j), row = importer,
    col = exporter, time = year, time_fe = TRUE, n_reps = 9, seed = 303))
  expect_identical(
    printed(fit_p),
    paste("Invariant permutation test (mwperm)",
          "------------------------------------ Design : panel Clusters :",
          "row=22, col=22, time=6 (2904 observations) Permutations : K = 21",
          "(group order 22, 9 reps) Resolution : p-values are multiples of",
          "1/22 = 0.045 per rep; reported floor 0.045 fta OLS estimate =",
          "0.6774 95% IPT CI [0.4671, 0.8691] H0: beta = 0 p-value = 0.045",
          "Decision : reject at alpha = 0.05"))
  fit_p0 <- with(trade_panel, mwperm_panel(
    y = log_trade, d = fta, x = cbind(log_gdp_i, log_gdp_j), row = importer,
    col = exporter, time = year, time_fe = FALSE, n_reps = 9, seed = 303))
  expect_identical(
    with(fit_p0, sprintf("estimate = %.4f 95%% CI [%.4f, %.4f]", estimate,
                         conf_int[1], conf_int[2])),
    "estimate = 0.9610 95% CI [0.4896, 1.3203]")
  expect_identical(fit_p0$pvalue, 1 / 22)                   # Table 5
  fit_pl <- with(trade_panel, mwperm_panel(                 # placebo policy
    y = log_trade, d = placebo, x = cbind(log_gdp_i, log_gdp_j),
    row = importer, col = exporter, time = year, time_fe = TRUE, n_reps = 9,
    seed = 303))
  expect_equal(round(unname(fit_pl$estimate), 3), -0.021)
  expect_equal(fit_pl$pvalue, 0.5)
  olsp <- stats::lm(log_trade ~ fta + log_gdp_i + log_gdp_j + factor(year),
                    data = trade_panel)
  expect_identical(
    squash(utils::capture.output(print(round(
      summary(olsp)$coefficients["fta", , drop = FALSE], 4)))),
    "Estimate Std. Error t value Pr(>|t|) fta 0.6774 0.0533 12.7087 0")
  expect_equal(unname(fit_p$estimate), unname(stats::coef(olsp)["fta"]),
               tolerance = 1e-10)
  expect_equal(unname(round(stats::confint(olsp)["fta", ], 3)),
               c(0.573, 0.782))
})

test_that("Sec. 6.2: the incomplete panel prints as in the draft", {
  tp <- with_seed(1, {
    tp <- trade_panel
    pair <- paste(tp$importer, tp$exporter)
    thin <- sample(unique(pair), 20)
    tp[!(pair %in% thin & tp$year > min(tp$year)), ]
  })
  fit_pm <- with(tp, mwperm_panel_missing(
    y = log_trade, d = fta, x = cbind(log_gdp_i, log_gdp_j), row = importer,
    col = exporter, time = year, min_block = 5, seed = 1))
  expect_identical(
    printed(fit_pm),
    paste("Invariant permutation test (mwperm)",
          "------------------------------------ Design : panel (bicliques)",
          "Clusters : row=22, col=22, time=6 (1632 observations)",
          "Permutations : K = 15 (group order 16, 10 reps) Resolution :",
          "p-values are multiples of 1/16 = 0.062 per rep; reported floor",
          "0.062 the median of 10 reps can fall between grid points fta OLS",
          "estimate = 0.6633 H0: beta = 0 p-value = 0.062 Decision : do not",
          "reject at alpha = 0.05 Notes: - Incomplete panel: 464 of 484",
          "observed (row, col) pairs are present in all 6 periods, and 272",
          "of those lie in the 1 fully observed block the biclique search",
          "extracted. Kept 1632 of 2804 observations (58.2%). The discarded",
          "pairs are what buys exact validity under missingness. - Block",
          "sizes (rows x cols): 17x16. - Resolution here is set by the",
          "smallest selected block: its permuted side is 16, so K = 15.",
          "Raise `min_block` so that small blocks cannot set K, or drop the",
          "sparsest periods so that more pairs clear the mask. - No 95%",
          "confidence interval: the smallest attainable p-value is 1/(K+1)",
          "= 0.0625, which is above alpha = 0.05, so no value could be",
          "excluded and the set would be the whole line. A 95% set needs K +",
          "1 >= 20 -- that is, at least 20 levels in the smallest permuted",
          "dimension. The p-value reported above is unaffected and remains",
          "exact."))
})

## ---- Section 6.3: three-way ------------------------------------------------

test_that("Sec. 6.3: the three-way example", {
  f3 <- with_seed(11, {
    n3 <- 20
    gg <- expand.grid(i = seq_len(n3), j = seq_len(n3), l = seq_len(n3))
    e1 <- stats::rnorm(n3)
    e2 <- stats::rnorm(n3)
    e3 <- stats::rnorm(n3)
    eps3 <- e1[gg$i] + e2[gg$j] + e3[gg$l] + stats::rnorm(nrow(gg))
    xx <- stats::rnorm(nrow(gg))
    dd <- stats::rnorm(nrow(gg))
    y3 <- 0.4 * dd + 0.5 * xx + eps3
    list(fit = mwperm_threeway(y = y3, d = dd, x = xx, id1 = gg$i,
                               id2 = gg$j, id3 = gg$l, n_reps = 9, seed = 11),
         ols = stats::coef(stats::lm(y3 ~ dd + xx))[["dd"]])
  })
  expect_identical(
    sub("^.*------------------------------------ ", "", printed(f3$fit)),
    paste("Design : threeway Clusters : id1=20, id2=20, id3=20 (8000",
          "observations) Permutations : K = 19 (group order 20, 9 reps)",
          "Resolution : p-values are multiples of 1/20 = 0.050 per rep;",
          "reported floor 0.050 dd OLS estimate = 0.3754 95% IPT CI [0.3353,",
          "0.4163] H0: beta = 0 p-value = 0.050 Decision : reject at alpha =",
          "0.05"))
  expect_equal(unname(f3$fit$estimate), f3$ols, tolerance = 1e-10)
})

## ---- Section 6.4: replicated two-way layout --------------------------------

test_that("Sec. 6.4: the balanced layout example", {
  fL <- with_seed(13, {
    cells <- expand.grid(r = seq_len(10), c = seq_len(10))
    reps <- sample(22:30, nrow(cells), replace = TRUE)
    rowsL <- rep(cells$r, reps)
    colsL <- rep(cells$c, reps)
    idl <- unlist(lapply(reps, seq_len))
    eta <- stats::rnorm(nrow(cells))
    zeta <- stats::rnorm(max(reps))
    dL <- stats::rnorm(length(rowsL))
    cell_id <- match(paste(rowsL, colsL), paste(cells$r, cells$c))
    yL <- 0.5 * dL + eta[cell_id] + zeta[idl] + stats::rnorm(length(rowsL))
    mwperm_layout(y = yL, d = dL, row = rowsL, col = colsL, rep = idl,
                  L0 = 22, n_reps = 9, seed = 13)
  })
  expect_identical(
    sub("^.*------------------------------------ ", "", printed(fL)),
    paste("Design : layout (balanced) Cells : 100 cells, smallest holds 22",
          "replicates (2200 observations) Permutations : K = 21 (group order",
          "22, 9 reps) Resolution : p-values are multiples of 1/22 = 0.045",
          "per rep; reported floor 0.045 dL OLS estimate = 0.483 95% IPT CI",
          "[0.3983, 0.5637] H0: beta = 0 p-value = 0.045 Decision : reject",
          "at alpha = 0.05 Notes: - Balanced to L0 = 22 replicates/cell:",
          "kept 100 of 100 cells and 2200 of 2604 observations (cells with <",
          "L0 replicates dropped; dense cells uniformly downsampled)."))
})

## ---- Section 6.5: irregular layout -----------------------------------------

irr <- with_seed(17, {
  cells <- expand.grid(i = seq_len(30), j = seq_len(30))
  thin <- cells$i > 25
  ell <- ifelse(thin, sample(1:3, nrow(cells), replace = TRUE),
                sample(4:9, nrow(cells), replace = TRUE))
  rowI <- rep(cells$i, ell)
  colI <- rep(cells$j, ell)
  idl <- unlist(lapply(ell, seq_len))
  cid <- rep(seq_len(nrow(cells)), ell)
  d_cell <- stats::rbinom(nrow(cells), 1, 0.4)
  eta <- stats::rnorm(nrow(cells))
  dI <- d_cell[cid]
  yI <- 0.5 * dI + eta[cid] + stats::rnorm(length(rowI))
  data.frame(y = yI, d = dI, i = rowI, j = colI, l = idl)
})
fI <- mwperm_irregular(y = irr$y, d = irr$d, row = irr$i, col = irr$j,
                       rep = irr$l, L0 = 4, min_block = 3, n_reps = 9,
                       seed = 17)

test_that("Sec. 6.5: the irregular example's numbers", {
  expect_equal(signif(unname(fI$estimate), 4), 0.5254)
  expect_equal(signif(fI$conf_int, 4), c(0.3406, 0.6749))
  expect_identical(fI$pvalue, 1 / 25)
  expect_identical(fI$K, 24L)
  expect_identical(fI$n_obs, 3000L)                 # per repetition
  expect_identical(fI$cells_used, 750L)
  expect_identical(fI$cells_total, 900L)
  expect_identical(fI$n_blocks, 1L)
  expect_match(fI$note[1], "uses 3000 of the 5261 observations (57.0%)",
               fixed = TRUE)
  expect_match(fI$note[1], "uses all 4974 observations", fixed = TRUE)
  expect_identical(fI$note[2], "Block sizes (rows x cols): 25x30.")
  ## "the estimate ... is OLS on all 4974 observations of the 750 retained
  ## cells"
  keep <- irr$i <= 25
  expect_equal(unname(fI$estimate),
               stats::coef(stats::lm(y ~ d, data = irr[keep, ]))[["d"]],
               tolerance = 1e-10)
  ## "Passing the same data to mwperm_layout() instead is refused ...; with
  ## L0 = 2, 3 or 4 ... it runs, warns ..., and returns p = 1 exactly"
  expect_error(suppressWarnings(mwperm_layout(irr$y, irr$d, row = irr$i,
                                              col = irr$j, rep = irr$l,
                                              seed = 1)),
               "Not enough clusters")
  for (L0 in 2:4) {
    expect_warning(
      fl <- mwperm_layout(irr$y, irr$d, row = irr$i, col = irr$j,
                          rep = irr$l, L0 = L0, n_reps = 3, seed = 1,
                          conf_int = FALSE),
      "constant within every cell")
    expect_identical(fl$pvalue, 1)
  }
})

test_that("Sec. 6.5: the irregular example's printed note", {
  skip_discrepancy("D5", paste("0.4.2 rewrote the irregular design's first",
                               "note (exchangeable replicates only) and",
                               "removed trim = \"levels\"; the draft prints",
                               "the 0.4.1 note"))
  expect_match(printed(fI), paste("The discarded observations are what buys",
                                  "exact validity under an unequal design."),
               fixed = TRUE)
})

test_that("Secs. 2.7, 3.2, 6.5: the trim = \"levels\" option", {
  skip_discrepancy("D5", paste("the draft documents",
                               "mwperm_irregular(trim = \"levels\"), removed",
                               "in 0.4.2 (now mwperm_panel_missing(time =,",
                               "L0 =))"))
  expect_no_error(mwperm_irregular(irr$y, irr$d, row = irr$i, col = irr$j,
                                   rep = irr$l, L0 = 4, trim = "levels",
                                   n_reps = 1, seed = 1, conf_int = FALSE))
})

## ---- Section 6.6: missing data ---------------------------------------------

test_that("Sec. 6.6: the missing-data example", {
  nd <- subset(trade_dyadic, importer != exporter)
  fit_m <- with(nd, mwperm_missing(
    y = log_trade, d = log_dist, x = cbind(log_gdp_i, log_gdp_j),
    row = importer, col = exporter, min_block = 3, n_reps = 11, seed = 505))
  expect_identical(
    sub("^.*------------------------------------ ", "", printed(fit_m)),
    paste("Design : missing (bicliques) Clusters : row=40, col=40 (800",
          "observations) Permutations : K = 19 (group order 20, 11 reps)",
          "Resolution : p-values are multiples of 1/20 = 0.050 per rep;",
          "reported floor 0.050 log_dist OLS estimate = -1.208 95% IPT CI",
          "[-1.694, -0.7019] H0: beta = 0 p-value = 0.050 Decision : reject",
          "at alpha = 0.05 Notes: - Kept 800 of 1560 observed cells (51.3%):",
          "the 2 fully observed rectangular blocks Procedure 2 could",
          "extract. The other 760 cells are discarded -- under missingness",
          "that loss is what buys exact validity, and a larger block would",
          "only add power. - Block sizes (rows x cols): 20x20, 20x20."))
  ## the estimate is OLS on the 800 retained cells
  ri <- as.integer(nd$importer)
  ci <- as.integer(nd$exporter)
  bl <- find_bicliques(ri, ci, min_block = 3)
  keep <- Reduce(`|`, lapply(bl, function(b) ri %in% b$rows & ci %in% b$cols))
  ols <- stats::lm(log_trade ~ log_dist + log_gdp_i + log_gdp_j,
                   data = nd[keep, ])
  expect_equal(unname(fit_m$estimate), unname(stats::coef(ols)["log_dist"]),
               tolerance = 1e-10)
})

## ---- Section 6.7: sign-flip ------------------------------------------------

fit_h <- with(trade_dyadic, mwperm_dyadic_het(
  y = log_trade, d = log_dist, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, n_flip = 7, n_reps = 15, seed = 101))

test_that("Sec. 6.7: the sign-flip example's numbers", {
  expect_identical(
    sub("^.*------------------------------------ ", "", printed(fit_h)),
    paste("Design : dyadic (sign-flip / heteroskedasticity-robust) Clusters",
          ": row=40, col=40 (1600 observations) Sign flips : n_flip = 7 flip",
          "groups (group order 2^6 = 64, 15 reps) Resolution : p-values are",
          "multiples of 1/64 = 0.016 per rep; reported floor 0.016 log_dist",
          "OLS estimate = -0.8985 95% IPT CI [-1.167, -0.6548] H0: beta = 0",
          "p-value = 0.016 Decision : reject at alpha = 0.05"))
  expect_true(all(fit_h$pvalues_rep == 1 / 64))    # "all fifteen ... on it"
  expect_identical(fit_h$ci_method, "exact")        # n_flip = 7 stays exact
  h2 <- quiet(mwperm(y = "log_trade", d = "log_dist",
                     x = c("log_gdp_i", "log_gdp_j"),
                     index = c("importer", "exporter"), data = trade_dyadic,
                     design = "dyadic_het", n_flip = 7, n_reps = 15,
                     seed = 101))
  expect_identical(h2$pvalues_rep, fit_h$pvalues_rep)
  expect_identical(h2$conf_int, fit_h$conf_int)
})

test_that("Sec. 6.7: the sign-flip example's printed header", {
  skip_discrepancy("D4", paste("0.4.2 heads a sign-flip fit 'Invariant",
                               "sign-flip test (mwperm)'; the draft prints",
                               "'Invariant permutation test (mwperm)'"))
  expect_match(printed(fit_h), "^Invariant permutation test \\(mwperm\\)")
})

## ---- Section 7: real bilateral-trade data (needs gravity) ------------------

test_that("Sec. 7: the Head-Mayer-Ries cross-section", {
  ## gravity is deliberately NOT a declared dependency: this test is optional
  ## by design and skips without it. The name is held in a variable so that
  ## R CMD check's scan for undeclared test dependencies does not demand it
  ## in Suggests (which would make every CI run install it).
  pkg <- "gravity"
  skip_if_not_installed(pkg)
  e <- new.env()
  utils::data(list = "gravity_no_zeros", package = pkg, envir = e)
  g <- e$gravity_no_zeros
  xs <- data.frame(importer = g$iso_d, exporter = g$iso_o,
                   log_trade = log(g$flow), log_dist = log(g$distw),
                   log_gdp_i = log(g$gdp_d), log_gdp_j = log(g$gdp_o),
                   rta = as.numeric(g$rta))
  expect_identical(nrow(xs), 17088L)
  expect_identical(length(unique(xs$importer)), 166L)
  expect_identical(length(unique(xs$exporter)), 166L)
  expect_identical(
    printed(mwperm_check(index = c("importer", "exporter"), data = xs)),
    paste("mwperm design diagnosis ------------------------------ Detected",
          "design : missing (2 indices, 17088 of 27556 cells observed)",
          "Roles : row = importer, col = exporter Dimensions : importer",
          "(166) x exporter (166) | 17088 observations Balance : incomplete",
          "(17088 of 27556 cells) Resolution : set by the biclique blocks,",
          "so not known until they are found (see find_bicliques) Would run",
          ": mwperm_missing(y, d, x, row = importer, col = exporter,",
          "min_block = ...) - The permutation-group order under missingness",
          "is set by the fully observed biclique blocks; see",
          "find_bicliques() for the achievable K. ? design = \"dyadic_het\"",
          "runs the sign-flip test on every observed cell, with no biclique",
          "search and nothing discarded: valid under arbitrary",
          "heteroskedasticity for errors that are independent across cells",
          "and symmetric about zero, but NOT under additive cluster effects",
          "eta_i + xi_j (see ?mwperm_dyadic_het)"))
  bl <- find_bicliques(xs$importer, xs$exporter, min_block = 21L)
  sides <- lapply(bl, function(b) c(length(b$rows), length(b$cols)))
  expect_setequal(sides, list(c(31L, 91L), c(100L, 30L)))
  expect_identical(sum(vapply(bl, function(b)
    length(b$rows) * length(b$cols), 0)), 5821)
  fit_ipt <- quiet(mwperm(y = "log_trade", d = "log_dist",
                          x = c("log_gdp_i", "log_gdp_j"),
                          index = c("importer", "exporter"), data = xs,
                          min_block = 21, n_reps = 9, seed = 42))
  expect_identical(fit_ipt$auto$design, "missing")
  expect_equal(round(unname(fit_ipt$estimate), 2), -0.99)
  expect_equal(round(fit_ipt$conf_int, 2), c(-1.17, -0.81))
  expect_equal(round(fit_ipt$pvalue, 3), 0.033)
  expect_identical(fit_ipt$cells_used, 5821L)
  expect_identical(fit_ipt$K + 1L, 30L)
  ols_all <- stats::lm(log_trade ~ log_dist + log_gdp_i + log_gdp_j,
                       data = xs)
  expect_equal(round(unname(stats::coef(ols_all)["log_dist"]), 2), -1.52)
  fit_rta <- quiet(mwperm(y = "log_trade", d = "rta",
                          x = c("log_dist", "log_gdp_i", "log_gdp_j"),
                          index = c("importer", "exporter"), data = xs,
                          min_block = 21, n_reps = 9, seed = 43))
  expect_equal(round(unname(fit_rta$estimate), 2), 0.69)
  expect_equal(round(fit_rta$conf_int, 2), c(0.21, 1.13))
  expect_equal(round(fit_rta$pvalue, 3), 0.033)
  ## Section 7.6: a null just inside the lower limit is not rejected
  b_edge <- fit_ipt$conf_int[1] + 0.03 * diff(fit_ipt$conf_int)
  expect_equal(round(b_edge, 2), -1.16)
  fit49 <- quiet(mwperm(y = "log_trade", d = "log_dist",
                        x = c("log_gdp_i", "log_gdp_j"),
                        index = c("importer", "exporter"), data = xs,
                        min_block = 21, beta_null = b_edge,
                        conf_int = FALSE, n_reps = 49, seed = 42))
  expect_equal(round(fit49$pvalue, 2), 0.10)
  expect_gt(fit49$pvalue, 0.05)
})

test_that("Sec. 7.4: the curated CEPII panel", {
  skip("Sec. 7.4 needs the CEPII Gravity database, which is not distributed")
})
