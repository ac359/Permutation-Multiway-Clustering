## mwperm_formula() -- the formula interface, y ~ d | x.
##
## This front end has NO statistical content: it assembles y, d and x from
## formula algebra and forwards everything to mwperm(). So the test is an
## identity test -- the fitted object must match the data interface exactly,
## field for field, modulo `call`. Any difference here is an assembly bug: a
## term landing in the wrong matrix, a transformation evaluated in the wrong
## frame, or the time role lost on the way through.
##
## Cross-cutting contracts live elsewhere; see tests/README.md.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

same_but_call <- function(a, b) {
  if (!identical(sort(names(a)), sort(names(b)))) return(FALSE)
  for (f in setdiff(names(a), "call")) {
    x <- a[[f]]
    y <- b[[f]]
    if (f == "auto" && is.list(x) && is.list(y)) x$call <- y$call <- NULL
    if (!identical(x, y)) return(FALSE)
  }
  TRUE
}

data(trade_dyadic)
data(trade_panel)

## ---- 1. dyadic: y ~ d | x1 + x2 == the data interface ---------------------
f1 <- mwperm_formula(log_trade ~ log_dist | log_gdp_i + log_gdp_j,
                     data = trade_dyadic, index = c("importer", "exporter"),
                     n_reps = 3, seed = 1, verbose = FALSE)
g1 <- mwperm(y = "log_trade", d = "log_dist", x = c("log_gdp_i", "log_gdp_j"),
             index = c("importer", "exporter"), data = trade_dyadic,
             n_reps = 3, seed = 1, verbose = FALSE)
stopifnot(same_but_call(f1, g1))

## ---- 2. no-nuisance form and transformed terms ----------------------------
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

## ---- 3. joint d > 1 and the panel time role -------------------------------
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

## ---- 4. the accessors read the formula's own names, and validation --------
stopifnot(identical(coef(f1), setNames(as.numeric(f1$estimate), "log_dist")),
          identical(nobs(f1), f1$n_obs), nobs(f1) == 1600L)
expect_err(mwperm_formula(~log_dist, data = trade_dyadic,
                          index = c("importer", "exporter")),
           "two-sided")

## ---- 5. missing values are refused by name, never dropped silently -------
## model.matrix() applies getOption("na.action") and drops incomplete rows on
## its own, while the outcome keeps every row; the mismatch used to surface as
## "`x` must have the same number of rows as `y`", naming an argument the
## caller never passed. The package contract is that incomplete data is an
## error naming what is missing, and that a row is never dropped silently --
## a dropped row makes a complete array incomplete without anyone noticing.
td_na <- trade_dyadic
td_na$log_gdp_i[3L] <- NA
m_x <- msg_of(mwperm_formula(log_trade ~ log_dist | log_gdp_i + log_gdp_j,
                             data = td_na, index = c("importer", "exporter"),
                             n_reps = 1, seed = 1, verbose = FALSE))
stopifnot(!is.na(m_x), grepl("log_gdp_i", m_x, fixed = TRUE),
          grepl("missing", m_x, fixed = TRUE),
          !grepl("same number of rows", m_x, fixed = TRUE))
td_na <- trade_dyadic
td_na$log_dist[c(2L, 5L)] <- NA
m_d <- msg_of(mwperm_formula(log_trade ~ log_dist | log_gdp_i + log_gdp_j,
                             data = td_na, index = c("importer", "exporter"),
                             n_reps = 1, seed = 1, verbose = FALSE))
stopifnot(!is.na(m_d), grepl("log_dist", m_d, fixed = TRUE),
          grepl("2 row", m_d, fixed = TRUE))
td_na <- trade_dyadic
td_na$log_trade[7L] <- NA
m_y <- msg_of(mwperm_formula(log_trade ~ log_dist, data = td_na,
                             index = c("importer", "exporter"),
                             n_reps = 1, seed = 1, verbose = FALSE))
stopifnot(!is.na(m_y), grepl("log_trade", m_y, fixed = TRUE))
## the same data with the NA row removed by the caller runs, on 1599 rows
f_ok <- mwperm_formula(log_trade ~ log_dist, data = td_na[-7L, ],
                       index = c("importer", "exporter"), design = "missing",
                       n_reps = 1, seed = 1, conf_int = FALSE, verbose = FALSE)
stopifnot(nobs(f_ok) <= 1599L)

passed("test-formula.R")
