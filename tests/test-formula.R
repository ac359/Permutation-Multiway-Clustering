## Formula-interface identity tests (adopted from the audit's Phase-7
## prototype harness): mwperm_formula() must reproduce the data interface
## EXACTLY -- identical objects modulo `call` -- because it only assembles
## inputs; the statistical core is untouched. Base-R stopifnot style; fast.
library(mwperm)

same_but_call <- function(a, b) {
  if (!identical(sort(names(a)), sort(names(b)))) return(FALSE)
  for (f in setdiff(names(a), "call")) {
    x <- a[[f]]; y <- b[[f]]
    if (f == "auto" && is.list(x) && is.list(y)) x$call <- y$call <- NULL
    if (!identical(x, y)) return(FALSE)
  }
  TRUE
}

data(trade_dyadic)
data(trade_panel)

## ---- 1. dyadic: y ~ d | x1 + x2 == the data interface ------------------------
f1 <- mwperm_formula(log_trade ~ log_dist | log_gdp_i + log_gdp_j,
                     data = trade_dyadic, index = c("importer", "exporter"),
                     n_reps = 3, seed = 1, verbose = FALSE)
g1 <- mwperm(y = "log_trade", d = "log_dist", x = c("log_gdp_i", "log_gdp_j"),
             index = c("importer", "exporter"), data = trade_dyadic,
             n_reps = 3, seed = 1, verbose = FALSE)
stopifnot(same_but_call(f1, g1))

## ---- 2. no-nuisance form and transformed terms --------------------------------
f2 <- mwperm_formula(log_trade ~ log_dist, data = trade_dyadic,
                     index = c("importer", "exporter"),
                     n_reps = 2, seed = 4, conf_int = FALSE, verbose = FALSE)
g2 <- with(trade_dyadic,
           mwperm_dyadic(log_trade, cbind(log_dist = log_dist),
                         row = importer, col = exporter,
                         n_reps = 2, seed = 4, conf_int = FALSE))
stopifnot(identical(f2$pvalue, g2$pvalue),
          identical(f2$estimate, g2$estimate))
f3 <- mwperm_formula(log_trade ~ log_dist | log_gdp_i + I(log_gdp_j^2),
                     data = trade_dyadic, index = c("importer", "exporter"),
                     n_reps = 2, seed = 4, conf_int = FALSE, verbose = FALSE)
g3 <- with(trade_dyadic,
           mwperm_dyadic(log_trade, cbind(log_dist = log_dist),
                         x = cbind(log_gdp_i, log_gdp_j^2),
                         row = importer, col = exporter,
                         n_reps = 2, seed = 4, conf_int = FALSE))
stopifnot(identical(f3$pvalue, g3$pvalue),
          identical(f3$estimate, g3$estimate))

## ---- 3. joint d > 1 and the panel time role -----------------------------------
f4 <- mwperm_formula(log_trade ~ log_dist + border | log_gdp_i + log_gdp_j,
                     data = trade_dyadic, index = c("importer", "exporter"),
                     n_reps = 1, seed = 2, conf_int = FALSE, verbose = FALSE)
stopifnot(identical(f4$d_names, c("log_dist", "border")),
          length(f4$estimate) == 2L)
f5 <- mwperm_formula(log_trade ~ fta | log_gdp_i + log_gdp_j,
                     data = trade_panel, index = c("importer", "exporter"),
                     time = "year", n_reps = 1, seed = 1, conf_int = FALSE,
                     verbose = FALSE)
g5 <- with(trade_panel,
           mwperm_panel(log_trade, cbind(fta = fta),
                        x = cbind(log_gdp_i, log_gdp_j),
                        row = importer, col = exporter, time = year,
                        n_reps = 1, seed = 1, conf_int = FALSE))
stopifnot(identical(f5$pvalue, g5$pvalue),
          identical(f5$estimate, g5$estimate),
          identical(f5$auto$design, "panel"))

## ---- 4. accessors + validation -------------------------------------------------
stopifnot(identical(coef(f1), setNames(as.numeric(f1$estimate), "log_dist")),
          identical(nobs(f1), f1$n_obs), nobs(f1) == 1600L)
msg <- tryCatch({ mwperm_formula(~log_dist, data = trade_dyadic,
                                 index = c("importer", "exporter")); NA },
                error = function(e) conditionMessage(e))
stopifnot(!is.na(msg), grepl("two-sided", msg))

cat("test-formula.R: all assertions passed\n")
