## ---------------------------------------------------------------------------
## Golden baseline for the numerical-invariance contract.
##
## Runs every exported mwperm_* front end plus mwperm() on the shipped data at
## default arguments with seed = 1, plus one non-default configuration per
## design (non-zero beta_null, n_reps = 9, conf_int = TRUE), plus
## find_bicliques() on the diagonal-deleted dyadic design, and saves the full
## result objects.
##
##   Rscript tests/golden/make_baseline.R          # write baseline.rds
##   Rscript tests/golden/make_baseline.R --check  # compare, field by field
##
## Requires a FRESHLY INSTALLED package (`R CMD INSTALL .`): several fits reach
## internals that resolve against the installed copy, not the working tree.
## ---------------------------------------------------------------------------

library(mwperm)
args <- commandArgs(trailingOnly = TRUE)
check <- "--check" %in% args
## Locate the snapshot directory from wherever the script is driven: the
## package root (Rscript tests/golden/make_baseline.R), the copied tests
## directory that R CMD check runs in, or the covr equivalent, which is
## named <pkg>-tests rather than tests -- so match on the directory that
## exists, never on the working directory's name.
here <- if (dir.exists(file.path("tests", "golden"))) {
  file.path("tests", "golden")
} else if (dir.exists("golden")) {
  "golden"
} else {
  "."
}
## --against=<file> compares against a different snapshot; used to diff the
## current tree against the pre-0.3.0 reference (baseline-0.2.0.rds) and see
## exactly which numbers the 0.3.0 changes moved.
against <- grep("^--against=", args, value = TRUE)
rds <- file.path(here, if (length(against))
                          sub("^--against=", "", against[1L])
                        else "baseline.rds")

data(trade_dyadic, package = "mwperm")
data(trade_panel,  package = "mwperm")

td <- trade_dyadic
tp <- trade_panel

## Incomplete dyadic array: drop self-trade, then thin to 90% of the remainder.
set.seed(42)
inc <- td[td$importer != td$exporter, ]
inc <- inc[sample.int(nrow(inc), floor(0.9 * nrow(inc))), ]

## Replicated two-way layout: 6 x 6 cells, 6 within-cell replicates.
set.seed(7)
cells <- expand.grid(i = 1:6, j = 1:6)
lay <- do.call(rbind, lapply(seq_len(nrow(cells)), function(r)
  data.frame(i = cells$i[r], j = cells$j[r], l = 1:6,
             d = rnorm(6), eta = rnorm(1))))
lay$y <- 0.6 * lay$d + lay$eta + rnorm(nrow(lay))

## Irregular layout for the Section 6.4 design: unequal cell sizes and a
## covariate that is CONSTANT within each cell (the case Section 6.4 exists
## for). 8 x 8 cells, some cells thin or empty.
set.seed(11)
irr <- do.call(rbind, lapply(seq_len(8L), function(i)
  do.call(rbind, lapply(seq_len(8L), function(j) {
    L <- sample(c(0L, 2L, 4L, 6L, 9L), 1L)
    if (L == 0L) return(NULL)
    data.frame(i = i, j = j, l = seq_len(L),
               d = rnorm(1),                       # constant within cell
               eta = rnorm(1))
  }))))
irr$y <- 0.5 * irr$d + irr$eta + rnorm(nrow(irr))

## Incomplete panel for the blockwise-InvB design: 20 country pairs are
## observed in the first year only, so they fail the "present in every period"
## mask and are dropped whole.
set.seed(23)
pp <- paste(tp$importer, tp$exporter)
thin_pairs <- sample(unique(pp), 20L)
ipn <- tp[!(pp %in% thin_pairs & tp$year > min(tp$year)), ]

B <- list()
run <- function(nm, expr) {
  B[[nm]] <<- tryCatch(withCallingHandlers(expr,
                         warning = function(w) invokeRestart("muffleWarning")),
                       error = function(e) paste("ERROR:", conditionMessage(e)))
  invisible(NULL)
}

## ---- defaults, seed = 1 ---------------------------------------------------
run("dyadic_default", with(td, mwperm_dyadic(
  y = log_trade, d = log_dist, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, seed = 1)))
run("panel_default", with(tp, mwperm_panel(
  y = log_trade, d = fta, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, time = year, seed = 1)))
run("threeway_default", with(tp, mwperm_threeway(
  y = log_trade, d = fta, x = cbind(log_gdp_i, log_gdp_j),
  id1 = importer, id2 = exporter, id3 = year, seed = 1)))
run("layout_default", with(lay, mwperm_layout(
  y = y, d = d, row = i, col = j, rep = l, seed = 1)))
run("missing_default", with(inc, mwperm_missing(
  y = log_trade, d = log_dist, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, min_block = 3L, seed = 1)))
run("panel_missing_default", with(ipn, mwperm_panel_missing(
  y = log_trade, d = fta, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, time = year, min_block = 5L, seed = 1)))
run("irregular_default", with(irr, mwperm_irregular(
  y = y, d = d, row = i, col = j, rep = l, L0 = 4L, min_block = 2L, seed = 1)))
run("unified_dyadic", mwperm(
  y = "log_trade", d = "log_dist", x = c("log_gdp_i", "log_gdp_j"),
  index = c("importer", "exporter"), data = td, seed = 1, verbose = FALSE))
run("unified_panel", mwperm(
  y = "log_trade", d = "fta", x = c("log_gdp_i", "log_gdp_j"),
  index = c("importer", "exporter"), time = "year", data = tp,
  seed = 1, verbose = FALSE))
run("formula_dyadic", mwperm_formula(
  log_trade ~ log_dist | log_gdp_i + log_gdp_j, data = td,
  index = c("importer", "exporter"), seed = 1, verbose = FALSE))
run("readme_quickstart", mwperm(
  y = "log_trade", d = "log_dist", x = c("log_gdp_i", "log_gdp_j"),
  index = c("importer", "exporter"), data = td, n_reps = 15, seed = 1,
  verbose = FALSE))

## ---- one non-default configuration per design -----------------------------
## non-zero beta_null, n_reps = 9, conf_int = TRUE
run("dyadic_nondefault", with(td, mwperm_dyadic(
  y = log_trade, d = log_dist, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, beta_null = -1, n_reps = 9L,
  conf_int = TRUE, seed = 1)))
run("panel_nondefault", with(tp, mwperm_panel(
  y = log_trade, d = fta, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, time = year, beta_null = 0.5, n_reps = 9L,
  conf_int = TRUE, seed = 1)))
run("threeway_nondefault", with(tp, mwperm_threeway(
  y = log_trade, d = fta, x = cbind(log_gdp_i, log_gdp_j),
  id1 = importer, id2 = exporter, id3 = year, beta_null = 0.5, n_reps = 9L,
  conf_int = TRUE, seed = 1)))
run("layout_nondefault", with(lay, mwperm_layout(
  y = y, d = d, row = i, col = j, rep = l, beta_null = 0.6, n_reps = 9L,
  conf_int = TRUE, seed = 1)))
run("missing_nondefault", with(inc, mwperm_missing(
  y = log_trade, d = log_dist, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, min_block = 3L, beta_null = -1,
  n_reps = 9L, conf_int = TRUE, seed = 1)))
run("panel_missing_nondefault", with(ipn, mwperm_panel_missing(
  y = log_trade, d = fta, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, time = year, min_block = 5L,
  beta_null = 0.5, n_reps = 9L, conf_int = TRUE, time_fe = FALSE, seed = 1)))
## Irregular design with a covariate that VARIES within cells and cells whose
## level sets differ (rows 1-20 hold periods 1-2, rows 21-40 periods 2-3): the
## case in which the slot the permutation holds fixed must be the `rep` level.
## 20 retained rows give K = 19, so the confidence set is attainable.
set.seed(3)
irr2 <- do.call(rbind, lapply(seq_len(40L), function(i)
  do.call(rbind, lapply(seq_len(20L), function(j)
    data.frame(i = i, j = j, t = if (i <= 20L) 1:2 else 2:3)))))
irr2$d <- as.numeric(irr2$t >= matrix(sample(1:4, 800L, TRUE), 40L,
                                      20L)[cbind(irr2$i, irr2$j)])
irr2$y <- 0.4 * irr2$d + rnorm(40L)[irr2$i] + rnorm(20L)[irr2$j] +
  c(0, 0, 5)[irr2$t] + rnorm(nrow(irr2))
run("irregular_nondefault", with(irr2, mwperm_irregular(
  y = y, d = d, row = i, col = j, rep = t, L0 = 2L, min_block = 2L,
  beta_null = 0.4, n_reps = 9L, conf_int = TRUE, seed = 1)))
run("layout_L0", with(irr, mwperm_layout(
  y = y, d = d, row = i, col = j, rep = l, L0 = 4L, n_reps = 9L,
  conf_int = TRUE, seed = 1)))
run("dyadic_multi", with(td, mwperm_dyadic(
  y = log_trade, d = cbind(log_dist, border), x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, n_reps = 9L, conf_int = TRUE, seed = 1)))
run("dyadic_grid", with(td, mwperm_dyadic(
  y = log_trade, d = log_dist, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, n_reps = 9L, conf_int = TRUE,
  grid = seq(-2, 0, by = 0.02), seed = 1)))
run("dyadic_placebo", with(td, mwperm_dyadic(
  y = log_trade, d = placebo, x = cbind(log_gdp_i, log_gdp_j),
  row = importer, col = exporter, conf_int = FALSE, seed = 1)))

## ---- find_bicliques() on the diagonal-deleted design ----------------------
diag_del <- td[td$importer != td$exporter, ]
run("bicliques_diag_greedy", with(diag_del, find_bicliques(
  row = importer, col = exporter, min_block = 3L, method = "greedy")))
run("bicliques_inc_greedy", with(inc, find_bicliques(
  row = importer, col = exporter, min_block = 3L, method = "greedy")))
run("bicliques_inc_min2", with(inc, find_bicliques(
  row = importer, col = exporter, min_block = 2L, method = "greedy")))
run("permset", build_perm_set(n = 20, K = 5, seed = 1))

## ---- write or check -------------------------------------------------------
## Fields added by later work (conf_set, aggregate, ...) are not in the
## baseline; comparison is over the baseline's own fields, so a purely
## ADDITIVE field never registers as drift.
if (!check) {
  saveRDS(B, rds)
  cat(sprintf("wrote %s (%d entries)\n", rds, length(B)))
} else {
  old <- readRDS(rds)
  bad <- 0L
  ## Two classes of field, compared two ways.
  ##
  ## Everything the inference actually reports -- the p-value (a count over
  ## K + 1, so exactly representable), the per-rep p-values, K, the cluster and
  ## observation counts, the labels and notes -- is reproducible to the bit on
  ## every platform, and is compared with identical(). Drift there is a defect.
  ##
  ## The fields below come out of BLAS: estimate and se_naive from the OLS
  ## solve, conf_int / conf_set / conf_region / conf_box from cross products
  ## and the breakpoint roots (u_j - v_k) / (M_j - W_k). A different BLAS sums
  ## the same terms in a different order, so these agree to ~1e-15 relative but
  ## not bitwise -- across Linux, Windows and macOS, and even between two
  ## builds on the same architecture. identical() on them tests the linear
  ## algebra library, not this package. The tolerance is still seven orders
  ## tighter than any real change: the 0.3.0 exact-CI rewrite moved endpoints
  ## by ~7e-3.
  approx_fields <- c("estimate", "se_naive", "conf_int", "conf_set",
                     "conf_region", "conf_box")
  same <- function(f, x, y) {
    if (f %in% approx_fields && is.numeric(x) && is.numeric(y))
      isTRUE(all.equal(x, y, tolerance = 1e-8))
    else identical(x, y)
  }
  for (nm in names(old)) {
    if (!nm %in% names(B)) { cat("MISSING: ", nm, "\n"); bad <- bad + 1L; next }
    a <- old[[nm]]; b <- B[[nm]]
    if (identical(a, b)) next
    if (is.list(a) && is.list(b) && !is.null(names(a))) {
      diff <- names(a)[!vapply(names(a),
                function(f) same(f, a[[f]], b[[f]]), logical(1))]
      diff <- setdiff(diff, "call")
      if (!length(diff)) next
      cat(sprintf("DIFFERS: %-22s fields: %s\n", nm,
                  paste(diff, collapse = ", ")))
      for (f in diff) {
        show1 <- function(v)
          paste(utils::capture.output(str(v)), collapse = " | ")
        cat("   ", f, " old: ", show1(a[[f]]), "\n")
        cat("   ", f, " new: ", show1(b[[f]]), "\n")
      }
    } else {
      cat(sprintf("DIFFERS: %-22s (not a named list)\n", nm))
    }
    bad <- bad + 1L
  }
  cat(if (bad == 0L) "baseline reproduces\n" else
      sprintf("%d baseline entries differ\n", bad))
}
