## Regression tests for the print-label clarification (OLS estimate / IPT CI)
## and the input-validation & display fixes. Base-R stopifnot style; fast.
library(mwperm)

msg_of <- function(expr) tryCatch({expr; NA_character_},
                                  error = function(e) conditionMessage(e))

## ---- shared toy data -------------------------------------------------------
set.seed(1)
n <- 21L                               # >= 20 clusters/dim so a 95% CI exists
g <- expand.grid(i = seq_len(n), j = seq_len(n)); N <- nrow(g)
d1 <- rnorm(n)[g$i] + rnorm(N)
y1 <- rnorm(n)[g$i] + rnorm(n)[g$j] + 0.4 * d1 + rnorm(N)

## ---- 1. print labels: "OLS estimate" and "IPT CI" / "IPT region" -----------
fit <- mwperm_dyadic(y1, d1, row = g$i, col = g$j, seed = 5)
out <- paste(capture.output(print(fit)), collapse = "\n")
stopifnot(grepl("OLS estimate =", out, fixed = TRUE))
stopifnot(grepl("% IPT CI [", out, fixed = TRUE))

D2 <- cbind(d1 = d1, d2 = rnorm(n)[g$j] + rnorm(N))
y2 <- y1 - 0.2 * D2[, 2L]
fit2 <- mwperm_dyadic(y2, D2, row = g$i, col = g$j, seed = 6)
out2 <- paste(capture.output(print(fit2)), collapse = "\n")
stopifnot(grepl("% IPT region [", out2, fixed = TRUE))

## returned values carry provenance: summary columns are prefixed by source
## (ols_/ipt_); confint keeps the generic's percentile labels but records its
## method in an attribute
tab <- suppressWarnings(capture.output(s <- summary(fit)))
stopifnot(identical(names(s), c("term", "ols_estimate", "ols_se_naive",
                                "ipt_ci_low", "ipt_ci_high", "p_value")))
ci <- confint(fit)
stopifnot(identical(colnames(ci), c("2.5 %", "97.5 %")))
stopifnot(identical(attr(ci, "method"), "IPT (inverted permutation test)"))
## the joint-region branch of confint carries the same attribute
stopifnot(identical(attr(confint(fit2), "method"),
                    "IPT (inverted permutation test)"))

## ---- 2. vector beta_null: one coherent H0 line, and plot() must not error --
fitj <- mwperm_dyadic(y2, D2, row = g$i, col = g$j, seed = 7,
                      beta_null = c(0.4, -0.2), conf_int = FALSE)
h0 <- grep("H0: beta =", capture.output(print(fitj)), value = TRUE)
stopifnot(length(h0) == 1L, grepl("0.4, -0.2", h0, fixed = TRUE))
grDevices::pdf(NULL)
plot(fitj)                             # errored pre-fix (barplot names.arg mismatch)
grDevices::dev.off()

## ---- 3. NA cluster ids error early with the argument's name ----------------
rowNA <- g$i; rowNA[1L] <- NA
m <- msg_of(mwperm_dyadic(y1, d1, row = rowNA, col = g$j, seed = 1))
stopifnot(grepl("`row`", m, fixed = TRUE), grepl("missing values", m))

## ---- 4. non-finite / non-numeric y, d, x error early by name ---------------
yy <- y1; yy[3L] <- NA
m <- msg_of(mwperm_dyadic(yy, d1, row = g$i, col = g$j, seed = 1))
stopifnot(grepl("`y`", m, fixed = TRUE), grepl("non-finite", m))
dd <- d1; dd[5L] <- Inf
m <- msg_of(mwperm_dyadic(y1, dd, row = g$i, col = g$j, seed = 1))
stopifnot(grepl("`d`", m, fixed = TRUE))
xx <- cbind(z = rnorm(N)); xx[7L, 1L] <- NaN
m <- msg_of(mwperm_dyadic(y1, d1, x = xx, row = g$i, col = g$j, seed = 1))
stopifnot(grepl("`x`", m, fixed = TRUE))
m <- msg_of(mwperm_dyadic(y1, rep(letters[1:3], length.out = N),
                          row = g$i, col = g$j, seed = 1))
stopifnot(grepl("`d` must be numeric", m, fixed = TRUE))

## layout reaches its own early check (the within-cell diagnostic reads `d`)
gl <- expand.grid(i = 1:4, j = 1:4, l = 1:5); Nl <- nrow(gl)
dl <- rnorm(Nl); dl[2L] <- NA
m <- msg_of(mwperm_layout(rnorm(Nl), dl, row = gl$i, col = gl$j, seed = 1))
stopifnot(grepl("`d`", m, fixed = TRUE), grepl("non-finite", m))

## ---- 5. scalar-argument validation ------------------------------------------
m <- msg_of(mwperm_dyadic(y1, d1, row = g$i, col = g$j, alpha = 0))
stopifnot(grepl("`alpha`", m, fixed = TRUE))
m <- msg_of(mwperm_dyadic(y1, d1, row = g$i, col = g$j, alpha = c(0.05, 0.1)))
stopifnot(grepl("`alpha`", m, fixed = TRUE))
m <- msg_of(mwperm_dyadic(y1, d1, row = g$i, col = g$j, n_reps = 0))
stopifnot(grepl("`n_reps`", m, fixed = TRUE))
m <- msg_of(mwperm_dyadic(y1, d1, row = g$i, col = g$j, n_cores = NA))
stopifnot(grepl("`n_cores`", m, fixed = TRUE))
m <- msg_of(mwperm_dyadic(y1, d1, row = g$i, col = g$j,
                          grid = c(0, NA, 1)))
stopifnot(grepl("`grid`", m, fixed = TRUE))

## beta_null must have length 1 or d (was silently recycled/truncated before)
m <- msg_of(mwperm_dyadic(y1, d1, row = g$i, col = g$j, beta_null = c(0, 1)))
stopifnot(grepl("`beta_null`", m, fixed = TRUE))
m <- msg_of(mwperm_dyadic(y1, d1, row = g$i, col = g$j, beta_null = NA))
stopifnot(grepl("`beta_null`", m, fixed = TRUE))

## ---- 6. valid seeded results are unchanged by the validation layer ----------
fit_again <- mwperm_dyadic(y1, d1, row = g$i, col = g$j, seed = 5)
stopifnot(identical(fit$pvalue, fit_again$pvalue),
          identical(fit$conf_int, fit_again$conf_int),
          identical(fit$estimate, fit_again$estimate))

cat("test-fixes.R: all assertions passed\n")
