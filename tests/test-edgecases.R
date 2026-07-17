## Edge-case / adversarial-input regression tests (audit Phase 5, brief §5).
## Base-R stopifnot style; fast. Pins the VALIDATION CONTRACT: every bad input
## fails early with an error naming the offending user-facing argument
## (CLAUDE.md §5.7), degenerate-but-legal designs stay valid, and the S3
## methods honour their documented guarantees. The behaviours the audit found
## deficient (silent factor-y coercion F5.1, noise-driven degenerate-d p-value
## F5.2, K/seed validation style F5.3/F2.3, retained-cells-only validation
## F5.4, NA in layout `rep` F5.7, ignored confint(parm=) F5.8) are FIXED and
## asserted below; the dispatch findings (F5.5/F3.4) live in
## test-detection.R.
library(mwperm)

msg_of <- function(expr) tryCatch({ expr; NA_character_ },
                                  error = function(e) conditionMessage(e))
expect_err <- function(expr, pattern) {
  m <- msg_of(expr)
  stopifnot(!is.na(m), grepl(pattern, m, fixed = TRUE))
}

## ---- shared fixtures --------------------------------------------------------
set.seed(101)
g6 <- expand.grid(i = 1:6, j = 1:6); N6 <- nrow(g6)
d6 <- rnorm(6)[g6$i] + rnorm(N6)
y6 <- rnorm(6)[g6$i] + rnorm(6)[g6$j] + 0.4 * d6 + rnorm(N6)
x6 <- rnorm(N6)
mkNA <- function(v, i = 3L, val = NA) { v[i] <- val; v }

## ---- 1. NA / NaN / Inf / non-numeric error EARLY, naming the argument -------
expect_err(mwperm_dyadic(mkNA(y6), d6, row = g6$i, col = g6$j), "`y`")
expect_err(mwperm_dyadic(mkNA(y6, val = NaN), d6, row = g6$i, col = g6$j), "`y`")
expect_err(mwperm_dyadic(mkNA(y6, val = Inf), d6, row = g6$i, col = g6$j), "`y`")
expect_err(mwperm_dyadic(y6, mkNA(d6), row = g6$i, col = g6$j), "`d`")
expect_err(mwperm_dyadic(y6, d6, x = mkNA(x6), row = g6$i, col = g6$j), "`x`")
expect_err(mwperm_dyadic(y6, d6, row = mkNA(g6$i), col = g6$j), "`row`")
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = mkNA(g6$j)), "`col`")
expect_err(mwperm_dyadic(y6, factor(round(d6)), row = g6$i, col = g6$j), "`d`")
expect_err(mwperm_dyadic(y6, as.character(round(d6)), row = g6$i, col = g6$j), "`d`")
expect_err(mwperm_dyadic(y6, d6, x = data.frame(a = x6, f = factor(g6$i)),
                         row = g6$i, col = g6$j), "`x`")
## non-numeric character y coerces to NA (with R's own warning) then errors on `y`
expect_err(suppressWarnings(
  mwperm_dyadic(letters[1 + (seq_len(N6) %% 5)], d6, row = g6$i, col = g6$j)), "`y`")
## the other front ends share the same validation layer -- one NA case each
gp <- expand.grid(i = 1:4, j = 1:4, t = 1:3)
yp <- rnorm(nrow(gp)); dp <- rnorm(nrow(gp))
expect_err(mwperm_panel(yp, dp, row = gp$i, col = gp$j, time = mkNA(gp$t)), "`time`")
expect_err(mwperm_threeway(yp, dp, id1 = gp$i, id2 = gp$j, id3 = mkNA(gp$t)), "`id3`")
g8 <- expand.grid(i = 1:8, j = 1:8); g8 <- g8[g8$i != g8$j, ]
y8 <- rnorm(nrow(g8)); d8 <- rnorm(nrow(g8))
## validation runs on the FULL supplied data, before the biclique step
## discards cells (F5.4 fix): an NA errors whether its cell is retained...
bl8v <- find_bicliques(g8$i, g8$j, min_block = 3)
i_ret <- which(g8$i %in% bl8v[[1]]$rows & g8$j %in% bl8v[[1]]$cols)[1L]
expect_err(mwperm_missing(mkNA(y8, i_ret), d8, row = g8$i, col = g8$j,
                          min_block = 3), "`y`")
## ... or discarded
in_any <- rep(FALSE, nrow(g8))
for (b in bl8v) in_any <- in_any | (g8$i %in% b$rows & g8$j %in% b$cols)
if (any(!in_any))
  expect_err(mwperm_missing(mkNA(y8, which(!in_any)[1L]), d8, row = g8$i,
                            col = g8$j, min_block = 3), "`y`")
## when the smallest block caps K below the alpha resolution, the fit's note
## names the block as the cause (F3.2 fix); attainable alpha -> no such note
fm05 <- mwperm_missing(y8, d8, row = g8$i, col = g8$j, min_block = 3,
                       seed = 1, conf_int = FALSE)
stopifnot(1 / (fm05$K + 1) > 0.05,               # fixture really is capped
          any(grepl("smallest selected block", fm05$note)))
fm40 <- mwperm_missing(y8, d8, row = g8$i, col = g8$j, min_block = 3,
                       seed = 1, conf_int = FALSE, alpha = 0.4)
stopifnot(!any(grepl("smallest selected block", fm40$note)))
gl <- expand.grid(l = 1:5, i = 1:4, j = 1:4)
yl <- rnorm(nrow(gl)); dl <- rnorm(16)[(gl$i - 1L) * 4L + gl$j] + rnorm(nrow(gl))
expect_err(mwperm_layout(mkNA(yl), dl, row = gl$i, col = gl$j), "`y`")
## factor y is rejected in every front end, never coerced to its level codes
## (audit F5.1: as.numeric(factor) is silent data corruption)
yf6 <- factor(round(y6)); yfp <- factor(round(yp))
yf8 <- factor(round(y8)); yfl <- factor(round(yl))
expect_err(mwperm_dyadic(yf6, d6, row = g6$i, col = g6$j), "`y`")
expect_err(mwperm_panel(yfp, dp, row = gp$i, col = gp$j, time = gp$t), "`y`")
expect_err(mwperm_threeway(yfp, dp, id1 = gp$i, id2 = gp$j, id3 = gp$t), "`y`")
expect_err(mwperm_missing(yf8, d8, row = g8$i, col = g8$j, min_block = 3), "`y`")
expect_err(mwperm_layout(yfl, dl, row = gl$i, col = gl$j), "`y`")
## logicals are documented as allowed
fl <- mwperm_dyadic(y6 > 0, d6 > 0, row = g6$i, col = g6$j, seed = 1,
                    conf_int = FALSE)
stopifnot(inherits(fl, "mwperm"))

## ---- 2. structural validation ------------------------------------------------
expect_err(mwperm_dyadic(c(y6, 1), c(d6, 1), row = c(g6$i, 1), col = c(g6$j, 1)),
           "one observation per (row, col) cell")
expect_err(mwperm_missing(c(y8, 1), c(d8, 1), row = c(g8$i, 1), col = c(g8$j, 2)),
           "one observation per (row, col) cell")
expect_err(mwperm_dyadic(y6, d6, x = matrix(rnorm(N6 * 18), N6),
                         row = g6$i, col = g6$j), "N > 2p")
expect_err(mwperm_dyadic(y6, d6, row = g6$i[-1], col = g6$j[-1]), "same length")
## non-square complete arrays are supported; K follows the smaller dimension
f59 <- mwperm_dyadic(rnorm(54), rnorm(54), row = rep(1:6, 9),
                     col = rep(1:9, each = 6), seed = 1, conf_int = FALSE)
stopifnot(f59$K == 5L,
          identical(unname(f59$n_clusters), c(6L, 9L)))
## unused factor levels are dropped, not counted as clusters
fuf <- mwperm_dyadic(y6, d6, row = factor(g6$i, levels = 1:10), col = g6$j,
                     seed = 1, conf_int = FALSE)
stopifnot(fuf$n_clusters[["row"]] == 6L)
## character / non-consecutive-integer ids give the same dense coding
fch <- mwperm_dyadic(y6, d6, row = letters[g6$i], col = g6$j, seed = 1,
                     conf_int = FALSE)
fnc <- mwperm_dyadic(y6, d6, row = c(10L, 20L, 35L, 40L, 70L, 99L)[g6$i],
                     col = g6$j, seed = 1, conf_int = FALSE)
fpl <- mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 1, conf_int = FALSE)
stopifnot(identical(fch$pvalue, fpl$pvalue), identical(fnc$pvalue, fpl$pvalue))

## ---- 3. degenerate-but-legal inputs -------------------------------------------
## constant y: statistic identically zero, p = 1 (never spuriously small)
fcy <- mwperm_dyadic(rep(2, N6), d6, row = g6$i, col = g6$j, seed = 1,
                     conf_int = FALSE)
stopifnot(fcy$pvalue == 1)
## constant d / d collinear with x: OLS estimate honestly NA (rank deficient);
## the test warns that beta is unidentified and returns p = 1 EXACTLY (the
## exact-arithmetic answer: a residualized-away d has a_k = 0, so the
## minorization gives p = 1; float noise must never decide instead -- F5.2 fix)
warn_of <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
  })
  w
}
wcd <- warn_of(fcd <- mwperm_dyadic(y6, rep(1, N6), row = g6$i, col = g6$j,
                                    seed = 1, conf_int = FALSE))
stopifnot(is.na(fcd$estimate), fcd$pvalue == 1,
          any(grepl("`d`", wcd, fixed = TRUE)))
wcl <- warn_of(fcl <- mwperm_dyadic(y6, 2 * x6 + 1, x = x6, row = g6$i,
                                    col = g6$j, seed = 1, conf_int = FALSE))
stopifnot(is.na(fcl$estimate), fcl$pvalue == 1,
          any(grepl("`d`", wcl, fixed = TRUE)))

## ---- 4. scalar-argument validation --------------------------------------------
for (a in list(0, 1, -0.1, 2, NA_real_, c(0.05, 0.1), "0.05"))
  expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, alpha = a), "`alpha`")
for (r in list(0L, 1.5, NA_integer_, c(2L, 3L)))
  expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, n_reps = r), "`n_reps`")
## K validation names the argument, incl. NA / vector / fractional / 0 (F5.3 fix)
for (k in list(NA, NA_integer_, c(3, 4), 2.5, "3", 0L))
  expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, K = k), "`K`")
## seed: type validation and a loud, `seed`-named overflow error (F2.3 fix);
## the last safe seed for this design keeps working (frozen scheme)
for (s in list(c(1, 2), "1", NA_real_, Inf))
  expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = s), "`seed`")
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 2147484,
                         n_reps = 1, conf_int = FALSE), "`seed`")
## n_reps = 1 explicitly: the pin is on the frozen sub-seed scheme itself
## (rep r uses seed + r - 1, so more reps move the boundary down by n_reps-1)
fs <- mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 2147482,
                    n_reps = 1, conf_int = FALSE)
stopifnot(inherits(fs, "mwperm"))
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, beta_null = Inf),
           "`beta_null`")
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, beta_null = c(0, 1)),
           "`beta_null`")
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, grid = c(0, NA, 1)),
           "`grid`")
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, K = 6), "K + 1")
## K = 0 blames the argument, not the data (F5.3: it used to say "Not enough
## clusters" as if the design were at fault)
expect_err(mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, K = 0), "`K`")
expect_err(mwperm_layout(yl, dl, row = gl$i, col = gl$j, L0 = 1), "`L0`")
expect_err(mwperm_layout(yl, dl, row = gl$i, col = gl$j, L0 = 99), "L0 = 99")
## a length-d vector null is legal and runs
fv <- mwperm_dyadic(y6, cbind(a = d6, b = x6), row = g6$i, col = g6$j,
                    beta_null = c(0.4, 0), seed = 1, conf_int = FALSE)
stopifnot(length(fv$beta_null) == 2L)

## ---- 5. resolution rule ---------------------------------------------------------
## K = 5 -> p_min = 1/6 > .05: p-value valid and on the grid, CI refused with a note
f6 <- mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 1)   # alpha = .05
stopifnot(is.null(f6$conf_int),
          any(grepl("1/(K+1)", f6$note, fixed = TRUE)),
          isTRUE(all.equal(f6$pvalue %% (1 / 6), 0)) ||
            f6$pvalue %in% ((1:6) / 6))
## same fit at alpha = .2 (attainable): the interval IS computed
f6b <- mwperm_dyadic(y6, d6, row = g6$i, col = g6$j, seed = 1, alpha = 0.2)
stopifnot(length(f6b$conf_int) == 2L, !anyNA(f6b$conf_int))

## ---- 6. confint contract (H3) ----------------------------------------------------
set.seed(2)
g21 <- expand.grid(i = 1:21, j = 1:21); N21 <- nrow(g21)
dA <- rnorm(21)[g21$i] + rnorm(N21)
dB <- rnorm(21)[g21$j] + rnorm(N21)
y21 <- rnorm(21)[g21$i] + rnorm(21)[g21$j] + 0.4 * dA - 0.2 * dB + rnorm(N21)
fit1 <- mwperm_dyadic(y21, dA, row = g21$i, col = g21$j, seed = 3)
expect_err(confint(fit1, level = 0.90), "refitting")     # wrong level: error, never mislabel
ci <- confint(fit1, level = 0.95)                        # matching level: fine
stopifnot(identical(colnames(ci), c("2.5 %", "97.5 %")),
          isTRUE(all.equal(as.numeric(ci), as.numeric(fit1$conf_int))),
          identical(attr(ci, "method"), "IPT (inverted permutation test)"))
fitn <- mwperm_dyadic(y21, dA, row = g21$i, col = g21$j, seed = 3, conf_int = FALSE)
expect_err(confint(fitn), "No confidence set")

## ---- 7. d > 1: summary/confint carry the joint region for EVERY term (H12) -------
fit2 <- mwperm_dyadic(y21, cbind(dA = dA, dB = dB), row = g21$i, col = g21$j,
                      seed = 3)
stopifnot(is.null(fit2$conf_int), !is.null(fit2$conf_box),
          nrow(fit2$conf_region) >= 1L)
o <- capture.output(tab <- summary(fit2))
stopifnot(nrow(tab) == 2L,
          isTRUE(all.equal(unname(tab$ipt_ci_low),  unname(fit2$conf_box[1, ]))),
          isTRUE(all.equal(unname(tab$ipt_ci_high), unname(fit2$conf_box[2, ]))),
          all(tab$p_value == fit2$pvalue))
ci2 <- confint(fit2)
stopifnot(nrow(ci2) == 2L, identical(rownames(ci2), c("dA", "dB")),
          isTRUE(all.equal(as.numeric(ci2), as.numeric(t(fit2$conf_box)))))
## confint(parm=) subsets by name or index, keeping labels and the method
## attribute (F5.8 fix; it was previously accepted and ignored)
cip <- confint(fit2, parm = "dB")
stopifnot(nrow(cip) == 1L, identical(rownames(cip), "dB"),
          identical(cip[1L, ], ci2[2L, ]),
          identical(attr(cip, "method"), attr(ci2, "method")))
cin <- confint(fit2, parm = 1L)
stopifnot(identical(rownames(cin), "dA"), identical(cin[1L, ], ci2[1L, ]))
expect_err(confint(fit2, parm = "nope"), "`parm`")
expect_err(confint(fit2, parm = 5), "`parm`")
expect_err(mwperm_dyadic(y21, cbind(dA, dB), row = g21$i, col = g21$j,
                         grid = list(1:3), seed = 3), "one vector per coefficient")
## d = 2 without a region: summary degrades to NA limits, no error
o <- capture.output(tabn <- summary(
  mwperm_dyadic(y21, cbind(dA = dA, dB = dB), row = g21$i, col = g21$j,
                seed = 3, conf_int = FALSE)))
stopifnot(all(is.na(tabn$ipt_ci_low)), all(is.na(tabn$ipt_ci_high)))

## ---- 8. print / summary / plot run on every design (incl. n_reps = 1) ------------
g3 <- expand.grid(i = 1:4, j = 1:4, k = 1:4)
y3 <- rnorm(4)[g3$i] + rnorm(4)[g3$j] + rnorm(nrow(g3)); d3 <- rnorm(nrow(g3))
fits <- list(
  dyadic   = fpl,
  panel    = mwperm_panel(yp, dp, row = gp$i, col = gp$j, time = gp$t,
                          seed = 1, conf_int = FALSE),
  threeway = mwperm_threeway(y3, d3, id1 = g3$i, id2 = g3$j, id3 = g3$k,
                             seed = 1, conf_int = FALSE),
  layout   = mwperm_layout(yl, dl, row = gl$i, col = gl$j, seed = 1,
                           conf_int = FALSE),
  missing  = mwperm_missing(y8, d8, row = g8$i, col = g8$j, min_block = 3,
                            seed = 1, conf_int = FALSE)
)
for (f in fits) {
  out <- capture.output(print(f))
  stopifnot(any(grepl("OLS estimate", out)), any(grepl("p-value", out)))
  o <- capture.output(s <- summary(f))
  stopifnot(is.data.frame(s), nrow(s) == length(f$estimate))
  pdf(NULL)
  suppressMessages(plot(f))                       # n_reps = 1: dot-strip path
  dev.off()
}

## ---- 9. .cell_code exactness guard (2^53) -----------------------------------------
stopifnot(is.numeric(mwperm:::.cell_code(cbind(1:3, 1:3))))
expect_err(mwperm:::.cell_code(rbind(c(2^18, 2^18, 2^18))), "2^53")

## ---- 10. layout: explicit within-cell-unique `rep` makes the fit row-order invariant
set.seed(77)
repv <- gl$l
f0 <- mwperm_layout(yl, dl, row = gl$i, col = gl$j, rep = repv, seed = 3,
                    conf_int = FALSE)
prm <- sample(nrow(gl))
f1 <- mwperm_layout(yl[prm], dl[prm], row = gl$i[prm], col = gl$j[prm],
                    rep = repv[prm], seed = 3, conf_int = FALSE)
stopifnot(identical(f0$pvalue, f1$pvalue),
          abs(f0$estimate - f1$estimate) < 1e-10)
expect_err(mwperm_layout(yl, dl, row = gl$i, col = gl$j, rep = repv[-1]),
           "`rep`")
## NA in `rep` is rejected like every other identifier (F5.7 fix; it was
## previously silently ranked last within its cell)
expect_err(mwperm_layout(yl, dl, row = gl$i, col = gl$j,
                         rep = mkNA(as.numeric(repv))), "`rep`")

cat("test-edgecases.R: all assertions passed\n")
