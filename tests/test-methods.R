## The S3 methods and the figures: print, summary, confint, coef, nobs, plot
## and mwperm_save.
##
## These are the package's whole user interface once a fit exists, and their
## labels are a PUBLIC CONTRACT: "OLS estimate" and "IPT CI"/"IPT region" in
## print(), the provenance-prefixed summary() column names, and confint()'s
## "2.5 %"/"97.5 %" dimnames with its "method" attribute. Renaming any of them
## breaks user code, so they are pinned here.
##
## Two properties matter beyond the labels. Every figure draws only from fields
## stored on the fit -- the test is never re-run -- and every drawing path
## restores par() and returns the fit invisibly. And no method may error on a
## fit that legitimately carries no confidence set: it falls back with a
## message, because a coarse design is a normal outcome, not a failure.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

grDevices::pdf(NULL)                        # null device: nothing hits disk

## ---- fixtures -------------------------------------------------------------
set.seed(42)
n <- 21L                            # >= 20 clusters/dim so a 95% CI exists
g <- expand.grid(i = seq_len(n), j = seq_len(n))
N <- nrow(g)
d1 <- rnorm(n)[g$i] + rnorm(N)
y1 <- rnorm(n)[g$i] + rnorm(n)[g$j] + 0.4 * d1 + rnorm(N)
D2 <- cbind(a = d1, b = rnorm(n)[g$j] + rnorm(N))
y2 <- y1 - 0.2 * D2[, "b"]

fit1 <- mwperm_dyadic(y1, d1, row = g$i, col = g$j, seed = 3, n_reps = 7)
fit2 <- mwperm_dyadic(y2, D2, row = g$i, col = g$j, seed = 4,
                      grid = list(seq(-0.5, 1.2, length.out = 11),
                                  seq(-1.0, 0.6, length.out = 11)))
stopifnot(!is.null(fit1$conf_int), nrow(fit2$conf_region) >= 1L)

## ---- 1. print(): the output-label contract --------------------------------
out1 <- paste(capture.output(print(fit1)), collapse = "\n")
stopifnot(grepl("OLS estimate =", out1, fixed = TRUE),
          grepl("% IPT CI [", out1, fixed = TRUE),
          grepl("p-value", out1))
out2 <- paste(capture.output(print(fit2)), collapse = "\n")
stopifnot(grepl("% IPT region [", out2, fixed = TRUE))

## a vector null prints ONE coherent H0 line listing the whole vector
fitj <- mwperm_dyadic(y2, D2, row = g$i, col = g$j, seed = 7,
                      beta_null = c(0.4, -0.2), conf_int = FALSE)
h0 <- grep("H0: beta =", capture.output(print(fitj)), value = TRUE)
stopifnot(length(h0) == 1L, grepl("0.4, -0.2", h0, fixed = TRUE))
plot(fitj)                          # errored pre-fix (barplot names.arg bug)

## the p-value formatter: NA prints as NA, tiny values do not print as 0
fmt_p <- internal(".fmt_p")
stopifnot(fmt_p(NA) == "NA", fmt_p(1e-5) == "< 0.001")

## ---- 2. summary(): a tidy data frame with provenance-prefixed columns -----
o <- capture.output(s1 <- summary(fit1))
stopifnot(is.data.frame(s1), nrow(s1) == 1L,
          identical(names(s1), c("term", "ols_estimate", "ols_se_naive",
                                 "ipt_ci_low", "ipt_ci_high", "p_value")))
## with d > 1 every row carries the joint region's marginal extent, and the
## p-value is the single joint p-value, repeated
o <- capture.output(s2 <- summary(fit2))
stopifnot(nrow(s2) == 2L,
          isTRUE(all.equal(unname(s2$ipt_ci_low),
                           unname(fit2$conf_box[1, ]))),
          isTRUE(all.equal(unname(s2$ipt_ci_high),
                           unname(fit2$conf_box[2, ]))),
          all(s2$p_value == fit2$pvalue))
## no region stored: the limits degrade to NA rather than erroring
o <- capture.output(s2n <- summary(
  mwperm_dyadic(y2, D2, row = g$i, col = g$j, seed = 3, conf_int = FALSE)))
stopifnot(all(is.na(s2n$ipt_ci_low)), all(is.na(s2n$ipt_ci_high)))

## ---- 3. confint(): dimnames, method attribute, level, parm ----------------
## The level is fixed when the test is inverted, so asking for another one must
## be an error -- relabelling a 95% set as 90% would be silently wrong.
expect_err(confint(fit1, level = 0.90), "refitting")
ci1 <- confint(fit1, level = 0.95)
stopifnot(identical(colnames(ci1), c("2.5 %", "97.5 %")),
          isTRUE(all.equal(as.numeric(ci1), as.numeric(fit1$conf_int))),
          identical(attr(ci1, "method"), "IPT (inverted permutation test)"))
expect_err(confint(mwperm_dyadic(y1, d1, row = g$i, col = g$j, seed = 3,
                                 conf_int = FALSE)),
           "No confidence set")
ci2 <- confint(fit2)
stopifnot(nrow(ci2) == 2L, identical(rownames(ci2), c("a", "b")),
          isTRUE(all.equal(as.numeric(ci2), as.numeric(t(fit2$conf_box)))))
## parm subsets by name or position, keeping the labels and the attribute
## (it used to be accepted and silently ignored)
cip <- confint(fit2, parm = "b")
stopifnot(nrow(cip) == 1L, identical(rownames(cip), "b"),
          identical(cip[1L, ], ci2[2L, ]),
          identical(attr(cip, "method"), attr(ci2, "method")))
cin <- confint(fit2, parm = 1L)
stopifnot(identical(rownames(cin), "a"), identical(cin[1L, ], ci2[1L, ]))
expect_err(confint(fit2, parm = "nope"), "`parm`")
expect_err(confint(fit2, parm = 5), "`parm`")

## ---- 4. coef() and nobs() -------------------------------------------------
stopifnot(identical(coef(fit1), setNames(as.numeric(fit1$estimate), "d1")),
          identical(nobs(fit1), fit1$n_obs), nobs(fit1) == N)
stopifnot(identical(names(coef(fit2)), c("a", "b")))

## ---- 5. every design prints, summarises and plots -------------------------
## Including at n_reps = 1, where the stability figure has a single point and
## takes the dot-strip path.
set.seed(101)
gp <- expand.grid(i = 1:6, j = 1:6, t = 1:3)
Np <- nrow(gp)
dp <- rnorm(6)[gp$i] + rnorm(Np)
yp <- rnorm(6)[gp$i] + rnorm(6)[gp$j] + c(0, 1, 2)[gp$t] + 0.5 * dp + rnorm(Np)
gt <- expand.grid(a = 1:5, b = 1:5, c = 1:4)
Nt <- nrow(gt)
dt3 <- rnorm(Nt)
yt3 <- rnorm(5)[gt$a] + rnorm(5)[gt$b] + rnorm(4)[gt$c] + 0.3 * dt3 + rnorm(Nt)
gl <- expand.grid(i = 1:5, j = 1:5, l = 1:4)
Nl <- nrow(gl)
dl <- rnorm(Nl)
yl <- rnorm(25L)[(gl$i - 1L) * 5L + gl$j] + 0.3 * dl + rnorm(Nl)
gm <- expand.grid(i = 1:10, j = 1:10)
gm <- gm[gm$i != gm$j, ]
Nm <- nrow(gm)
dm <- rnorm(10)[gm$i] + rnorm(Nm)
ym <- rnorm(10)[gm$i] + rnorm(10)[gm$j] + 0.5 * dm + rnorm(Nm)

fits <- list(
  dyadic   = fit1,
  panel    = mwperm_panel(yp, dp, row = gp$i, col = gp$j, time = gp$t,
                          seed = 12, conf_int = FALSE),
  threeway = mwperm_threeway(yt3, dt3, id1 = gt$a, id2 = gt$b, id3 = gt$c,
                             seed = 13, conf_int = FALSE),
  layout   = mwperm_layout(yl, dl, row = gl$i, col = gl$j, seed = 14,
                           conf_int = FALSE),
  missing  = mwperm_missing(ym, dm, row = gm$i, col = gm$j, min_block = 3,
                            seed = 11, conf_int = FALSE)
)
for (f in fits) {
  out <- capture.output(print(f))
  stopifnot(any(grepl("OLS estimate", out)), any(grepl("p-value", out)))
  o <- capture.output(sm <- summary(f))
  stopifnot(is.data.frame(sm), nrow(sm) == length(f$estimate))
  vis <- withVisible(suppressMessages(plot(f)))
  stopifnot(!vis$visible, identical(vis$value, f))  # returns the fit invisibly
  suppressMessages(plot(f, type = "all"))
}

## ---- 6. plot(): every type draws, and par() is restored -------------------
mar0 <- graphics::par("mar")
mfrow0 <- graphics::par("mfrow")
for (tp in c("auto", "coef", "stability", "all")) plot(fit1, type = tp)
stopifnot(identical(graphics::par("mar"), mar0),
          identical(graphics::par("mfrow"), mfrow0))
for (tp in c("auto", "coef", "region", "stability", "all")) plot(fit2,
                                                                 type = tp)
## a vector null draws per-row null crosses instead of one shared line
plot(mwperm_dyadic(y2, D2, row = g$i, col = g$j, seed = 4,
                   beta_null = c(0.4, -0.2),
                   grid = list(seq(-0.5, 1.2, length.out = 11),
                               seq(-1.0, 0.6, length.out = 11))),
     type = "coef")

## ---- 7. plot() falls back with a message, never errors --------------------
## An empty region (a grid that misses the acceptance set), d = 3 (no plane to
## draw), no stored confidence set, and the reserved type names.
fit_e <- mwperm_dyadic(y2, D2, row = g$i, col = g$j, seed = 8,
                       grid = list(c(100, 101), c(100, 101)))
plot(fit_e)
stopifnot(length(msgs_of(plot(fit_e, type = "region"))) >= 1L)
D3 <- cbind(D2, c = rnorm(N))
fit3 <- mwperm_dyadic(y2 + 0.1 * D3[, "c"], D3, row = g$i, col = g$j, seed = 5,
                      grid = list(seq(0.0, 0.8, length.out = 5),
                                  seq(-0.6, 0.2, length.out = 5),
                                  seq(-0.3, 0.5, length.out = 5)))
plot(fit3)
stopifnot(any(grepl("two coefficients", msgs_of(plot(fit3, type = "region")))))
fit_no <- mwperm_dyadic(y1, d1, row = g$i, col = g$j, seed = 6,
                        conf_int = FALSE)
stopifnot(any(grepl("stability", msgs_of(plot(fit_no, type = "coef")))))
plot(fit_no)                                # auto -> stability, silently fine
g6 <- expand.grid(i = 1:6, j = 1:6)         # coarse: conf_int NULL by design
fit6 <- mwperm_dyadic(rnorm(36), rnorm(36), row = g6$i, col = g6$j, seed = 7)
for (tp in c("auto", "coef", "stability", "all")) plot(fit6, type = tp)
stopifnot(length(msgs_of(plot(fit1, type = "null"))) >= 1L,
          length(msgs_of(plot(fit1, type = "profile"))) >= 1L)
## few reps: the stability figure switches to the dot strip
plot(mwperm_dyadic(y1, d1, row = g$i, col = g$j, seed = 1, n_reps = 2,
                   conf_int = FALSE), type = "stability")

## ---- 8. the style layer: Okabe-Ito, never colour alone, overridable -------
sty <- internal(".mwperm_style")
s <- sty()
stopifnot(identical(s$col_estimate, "#0072B2"),      # Okabe-Ito blue
          identical(s$col_null, "#D55E00"),          # Okabe-Ito vermillion
          s$pch_estimate != s$pch_null,              # marks differ, not just
          s$lty_null != s$lty_alpha)                 #   colour
s2 <- sty(col_estimate = "black", lwd_interval = 3)
stopifnot(identical(s2$col_estimate, "black"), s2$lwd_interval == 3)
stopifnot(grepl("no_such_element", msg_of(sty(no_such_element = 1)),
                fixed = TRUE))
plot(fit1, col_estimate = "black", lwd_interval = 3)
plot(fit1, main = "custom title", sub = "", xlab = "beta")
expect_warn(plot(fit1, bogus_argument = 1), "bogus_argument")

## ---- 9. mwperm_save(): journal dimensions on every device -----------------
for (ext in c("pdf", "png", "tiff", "jpeg")) {
  f <- file.path(tempdir(), paste0("mwperm-test-fig.", ext))
  out <- mwperm_save(fit1, f)
  stopifnot(identical(out, f), file.exists(f), file.size(f) > 0)
  unlink(f)
}
f_all <- file.path(tempdir(), "mwperm-test-all.pdf")
mwperm_save(fit1, f_all, width = "double", type = "all", pointsize = 9)
stopifnot(file.exists(f_all))
unlink(f_all)
expect_err(mwperm_save(fit1, 1), "`file`")
expect_err(mwperm_save(fit1, file.path(tempdir(), "f.bmp")), "extension")
expect_err(mwperm_save(fit1, "a.png", width = -1), "`width`")
expect_err(mwperm_save(fit1, "a.png", height = 0), "`height`")
expect_err(mwperm_save(fit1, "a.png", res = 10), "`res`")
expect_err(mwperm_save(list(), file.path(tempdir(), "f.png")), "mwperm")

invisible(grDevices::dev.off())
passed("test-methods.R")
