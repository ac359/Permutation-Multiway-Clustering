## The input-validation contract, shared by every front end.
##
## Validation is centralised (.check_finite / .dense_id / .check_lengths in
## R/engine.R, plus each front end's own early checks), so it is tested once
## here rather than five times. The contract has three parts:
##
##   * bad input fails EARLY, with an error naming the offending user-facing
##     argument in backticks -- never a coercion error from deep inside
##     lm.fit() or a QR, and never a silent coercion (a factor `y` turned into
##     its level codes is data corruption that no downstream check can catch);
##   * degenerate-but-legal input still produces the exact-arithmetic answer,
##     not one decided by floating-point noise;
##   * when the design cannot support the requested alpha, the package refuses
##     to report a confidence set and says why, instead of returning one that
##     looks informative.
##
## Argument checks that belong to a single design -- `L0`, `rep`, `permute`,
## the one-observation-per-cell rules -- are in that design's own file. The
## data-side checks (lengths, `N > 2p`, the vector, column-name and formula
## interfaces, cluster-id types) are in tests/testthat/test-data-validation.R.
## See tests/README.md.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

## ---- fixtures: one small design of each shape -----------------------------
set.seed(101)
g6 <- expand.grid(i = 1:6, j = 1:6)                 # dyadic
N6 <- nrow(g6)
d6 <- rnorm(6)[g6$i] + rnorm(N6)
y6 <- rnorm(6)[g6$i] + rnorm(6)[g6$j] + 0.4 * d6 + rnorm(N6)
x6 <- rnorm(N6)

gp <- expand.grid(i = 1:4, j = 1:4, t = 1:3)        # panel / threeway
yp <- rnorm(nrow(gp))
dp <- rnorm(nrow(gp))

g8 <- expand.grid(i = 1:8, j = 1:8)                 # incomplete
g8 <- g8[g8$i != g8$j, ]
y8 <- rnorm(nrow(g8))
d8 <- rnorm(nrow(g8))

gl <- expand.grid(l = 1:5, i = 1:4, j = 1:4)        # replicated layout
yl <- rnorm(nrow(gl))
dl <- rnorm(16)[(gl$i - 1L) * 4L + gl$j] + rnorm(nrow(gl))

mkNA <- function(v, i = 3L, val = NA) {
  v[i] <- val
  v
}

## ---- 1. non-finite and non-numeric input, named by argument ---------------
expect_err(mwperm_dyadic(mkNA(y6), d6, row = g6$i, col = g6$j), "`y`")
expect_err(mwperm_dyadic(mkNA(y6, val = NaN), d6, row = g6$i, col = g6$j),
           "`y`")
expect_err(mwperm_dyadic(mkNA(y6, val = Inf), d6, row = g6$i, col = g6$j),
           "`y`")
expect_err(mwperm_dyadic(y6, mkNA(d6), row = g6$i, col = g6$j), "`d`")
expect_err(mwperm_dyadic(y6, d6, x = mkNA(x6), row = g6$i, col = g6$j), "`x`")
expect_err(mwperm_dyadic(y6, d6, row = mkNA(g6$i), col = g6$j), "`row`")
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = mkNA(g6$j)), "`col`")
expect_err(mwperm_dyadic(y6, factor(round(d6)), row = g6$i, col = g6$j), "`d`")
expect_err(mwperm_dyadic(y6, as.character(round(d6)), row = g6$i, col = g6$j),
           "`d`")
expect_err(mwperm_dyadic(y6, d6, x = data.frame(a = x6, f = factor(g6$i)),
                         row = g6$i, col = g6$j), "`x`")
## a character y coerces to NA (with R's own warning) and then errors on `y`
expect_err(suppressWarnings(
  mwperm_dyadic(letters[1 + (seq_len(N6) %% 5)], d6, row = g6$i, col = g6$j)),
  "`y`")
## the message says WHAT is wrong, not just which argument
stopifnot(grepl("non-finite",
                msg_of(mwperm_dyadic(mkNA(y6), d6, row = g6$i, col = g6$j))),
          grepl("missing values",
                msg_of(mwperm_dyadic(y6, d6, row = mkNA(g6$i), col = g6$j))))

## every front end shares the layer: one non-finite case each, on its own
## design-specific identifier
expect_err(mwperm_panel(yp, dp, row = gp$i, col = gp$j, time = mkNA(gp$t)),
           "`time`")
expect_err(mwperm_threeway(yp, dp, id1 = gp$i, id2 = gp$j, id3 = mkNA(gp$t)),
           "`id3`")
expect_err(mwperm_missing(mkNA(y8), d8, row = g8$i, col = g8$j, min_block = 3),
           "`y`")
expect_err(mwperm_layout(mkNA(yl), dl, row = gl$i, col = gl$j), "`y`")
## the layout front end reads `d` in its own early diagnostic, before the
## engine's central check -- it must reject there too, not read the bad values
stopifnot(grepl("non-finite",
                msg_of(mwperm_layout(yl, mkNA(dl), row = gl$i, col = gl$j))))

## a factor `y` is rejected everywhere, never coerced to its level codes
expect_err(mwperm_dyadic(factor(round(y6)), d6, row = g6$i, col = g6$j), "`y`")
expect_err(mwperm_panel(factor(round(yp)), dp, row = gp$i, col = gp$j,
                        time = gp$t), "`y`")
expect_err(mwperm_threeway(factor(round(yp)), dp, id1 = gp$i, id2 = gp$j,
                           id3 = gp$t), "`y`")
expect_err(mwperm_missing(factor(round(y8)), d8, row = g8$i, col = g8$j,
                          min_block = 3), "`y`")
expect_err(mwperm_layout(factor(round(yl)), dl, row = gl$i, col = gl$j), "`y`")
## logicals ARE documented as allowed
stopifnot(inherits(mwperm_dyadic(y6 > 0, d6 > 0, row = g6$i, col = g6$j,
                                 seed = 1, conf_int = FALSE), "mwperm"))

## ---- 2. structural validation ---------------------------------------------
expect_err(mwperm_dyadic(y6, d6, row = g6$i[-1], col = g6$j[-1]), "same length")
## the stacked projection consumes up to 2p columns, so N > 2p is required
expect_err(mwperm_dyadic(y6, d6, x = matrix(rnorm(N6 * 18), N6), row = g6$i,
                         col = g6$j), "N > 2p")

## ---- 3. scalar arguments blame the argument, not the data -----------------
for (a in list(0, 1, -0.1, 2, NA_real_, c(0.05, 0.1), "0.05"))
  expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, alpha = a),
             "`alpha`")
for (r in list(0L, 1.5, NA_integer_, c(2L, 3L)))
  expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, n_reps = r),
             "`n_reps`")
## K = 0 used to report "Not enough clusters", as if the design were at fault
for (k in list(NA, NA_integer_, c(3, 4), 2.5, "3", 0L))
  expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, K = k), "`K`")
for (s in list(c(1, 2), "1", NA_real_, Inf))
  expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = s), "`seed`")
## the seed scheme is rep r -> seed + r - 1, then .sub_seed(., j): a seed near
## the integer ceiling overflows, and must say `seed` rather than fail obscurely
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 2147484,
                         n_reps = 1, conf_int = FALSE), "`seed`")
stopifnot(inherits(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j,
                                 seed = 2147482, n_reps = 1,
                                 conf_int = FALSE), "mwperm"))
## beta_null must have length 1 or d (it was silently recycled before)
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, beta_null = Inf),
           "`beta_null`")
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, beta_null = NA),
           "`beta_null`")
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, beta_null = c(0, 1)),
           "`beta_null`")
stopifnot(length(mwperm_dyadic(y6, cbind(a = d6, b = x6), row = g6$i,
                               col = g6$j, beta_null = c(0.4, 0), seed = 1,
                               conf_int = FALSE)$beta_null) == 2L)
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, grid = c(0, NA, 1)),
           "`grid`")
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, n_cores = NA),
           "`n_cores`")
expect_err(mwperm_dyadic(y6, cbind(d6, x6), row = g6$i, col = g6$j,
                         grid = list(1:3), seed = 3, alpha = 0.4),
           "one vector per coefficient")

## ---- 4. degenerate but legal: the exact answer, and an honest estimate ----
## A constant outcome makes the statistic identically zero. p must be exactly
## 1 -- never a small number produced by float noise in a zero-vs-zero
## comparison.
stopifnot(mwperm_dyadic(rep(2, N6), d6, row = g6$i, col = g6$j, seed = 1,
                        conf_int = FALSE)$pvalue == 1)
## d constant, or exactly collinear with x: beta is unidentified, so the OLS
## estimate is honestly NA and the test warns. A residualized-away d has
## a_k = 0 for every k, so the minorization gives p = 1 in exact arithmetic.
w_c <- warns_of(f_c <- mwperm_dyadic(y6, rep(1, N6), row = g6$i, col = g6$j,
                                     seed = 1, conf_int = FALSE))
stopifnot(is.na(f_c$estimate), f_c$pvalue == 1,
          any(grepl("`d`", w_c, fixed = TRUE)))
w_l <- warns_of(f_l <- mwperm_dyadic(y6, 2 * x6 + 1, x = x6, row = g6$i,
                                     col = g6$j, seed = 1, conf_int = FALSE))
stopifnot(is.na(f_l$estimate), f_l$pvalue == 1,
          any(grepl("`d`", w_l, fixed = TRUE)))
## with conf_int = TRUE the interval search has nothing to bracket: it must
## return the whole line rather than an arbitrary finite-looking interval
w_d <- warns_of(f_d <- mwperm_dyadic(rnorm(N6), rep(1, N6), row = g6$i,
                                     col = g6$j, seed = 1, n_reps = 1,
                                     alpha = 0.4))
stopifnot(any(grepl("`d`", w_d, fixed = TRUE)), f_d$pvalue == 1,
          all(is.infinite(f_d$conf_int)))

## ---- 5. the resolution rule: no confidence set it cannot support ----------
## p-values live on the grid {1..K+1}/(K+1), so a 95% set needs 1/(K+1) <= .05,
## i.e. at least 20 levels in the smallest permuted dimension. Below that the
## p-value is still exact; the interval is withheld with an explanation.
f_lo <- mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 1)  # K = 5
stopifnot(is.null(f_lo$conf_int),
          any(grepl("1/(K+1)", f_lo$note, fixed = TRUE)),
          f_lo$pvalue %in% ((1:6) / 6))
expect_err(confint(f_lo), "No confidence set")
## the same fit at an attainable alpha does compute the interval
f_hi <- mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 1, alpha = 0.2)
stopifnot(length(f_hi$conf_int) == 2L, !anyNA(f_hi$conf_int))
## The gate is the smallest REPORTED p-value, not the grid step. `median2`
## reports min(1, 2 * median), so its floor is 2/(K+1): at alpha = 0.2 the
## default rule can just reject and gets an interval, median2 cannot and is
## refused one on exactly the same data.
f_m1 <- mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 1, alpha = 0.2)
f_m2 <- mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 1, alpha = 0.2,
                      aggregate = "median2")
stopifnot(identical(f_m1$resolution, f_m2$resolution),
          identical(f_m1$resolution, 1 / 6),
          identical(f_m1$p_floor, 1 / 6), identical(f_m2$p_floor, 2 / 6),
          !is.null(f_m1$conf_int), is.null(f_m2$conf_int),
          any(grepl("1/(K+1)", f_m2$note, fixed = TRUE)))

## ---- 6. an EMPTY acceptance set is reported as such, never as [NA, NA] ---
## A grid that misses the acceptance set entirely (here 50..60, far from an
## estimate near 0.4) retains nothing. The interval is then NA and the set has
## no components -- a meaningful, alarming outcome -- and the fit must say so
## in a note that names `grid`, rather than print "[NA, NA]" with no comment.
set.seed(5)
g21 <- expand.grid(i = 1:21, j = 1:21)
d21 <- rnorm(21)[g21$i] + rnorm(nrow(g21))
y21 <- rnorm(21)[g21$i] + rnorm(21)[g21$j] + 0.4 * d21 + rnorm(nrow(g21))
f_empty <- mwperm_dyadic(y21, d21, row = g21$i, col = g21$j, seed = 1,
                         n_reps = 2, grid = seq(50, 60, by = 0.5))
stopifnot(all(is.na(f_empty$conf_int)), nrow(f_empty$conf_set) == 0L,
          identical(f_empty$ci_method, "grid"),
          any(grepl("EMPTY", f_empty$note, fixed = TRUE) &
              grepl("`grid`", f_empty$note, fixed = TRUE)))
## and print() says so next to the interval instead of showing [NA, NA]
o_empty <- capture.output(print(f_empty))
stopifnot(any(grepl("IPT CI: empty set", o_empty, fixed = TRUE)),
          !any(grepl("[NA, NA]", o_empty, fixed = TRUE)))

## ---- 7. n_cores is validated like every other scalar argument ------------
## 0, a negative, and a fraction used to run silently (and serially): a typo
## passed unremarked, against the documented "single integer >= 1".
for (nc in list(0L, -1L, 1.5, c(1L, 2L), "2", Inf))
  expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, n_cores = nc,
                           conf_int = FALSE, seed = 1),
             "`n_cores`")
stopifnot(inherits(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, n_cores = 1,
                                 conf_int = FALSE, seed = 1), "mwperm"))

passed("test-validation.R")
