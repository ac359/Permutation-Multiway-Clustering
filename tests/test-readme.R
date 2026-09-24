## The output README.md SHOWS must be the output the package PRODUCES.
##
## README.md displays five transcripts -- the quick-start dyadic fit and its
## summary() data frame, the design diagnosis, and the panel and sign-flip fits
## under "Extensions" -- as fenced blocks a reader is invited to reproduce.
## Nothing used to check them, and two of the three the README had at 0.3.0
## silently went stale for a whole release: they still carried 0.2.0
## intervals after the 0.3.0 exact confidence set moved eight seeded end
## points. tests/golden/ pins the fitted OBJECTS, which is why it did not
## catch this; what rotted was the printed transcript.
##
## Section 1 pins the exact printed lines and always runs, so a change in the
## package's output fails here whatever the working directory. Section 2 then
## checks those same literals actually appear in README.md, which closes the
## gap the other way -- an edit to this file that is not mirrored in the README
## -- but only when README.md can be found. It cannot be under R CMD check: the
## checker copies tests/ into <pkg>.Rcheck/tests and leaves the source tree
## behind, and README.md is not installed. Run from the package root to get it:
##
##   R CMD INSTALL . && Rscript tests/test-readme.R
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

data(trade_dyadic)
data(trade_panel)

## The lines README.md shows. Keep these in step with the README; if the
## package's numbers move, BOTH have to be updated, and section 2 says so.
readme_lines <- c(
  quick_start_est =
    "  log_dist     OLS estimate = -0.8985   95% IPT CI [-1.246, -0.5485]",
  quick_start_p =
    "H0: beta = 0    p-value = 0.025",
  panel_est =
    "  fta          OLS estimate = 0.6774   95% IPT CI [0.442, 0.8803]",
  panel_p =
    "H0: beta = 0    p-value = 0.045",
  check_resolution =
    "Resolution      : default K = 39, so p-values are multiples of 1/40 = 0.025",
  check_verdict =
    "                  -> fine enough for a 95% confidence set at alpha = 0.05",
  quick_start_res =
    "Resolution   : p-values are multiples of 1/40 = 0.025 per rep; reported floor 0.025",
  panel_res =
    "Resolution   : p-values are multiples of 1/22 = 0.045 per rep; reported floor 0.045",
  panel_res_caveat =
    "               the median of 10 reps can fall between grid points",
  quick_start_summary =
    "1 log_dist   -0.8985026   0.08728893  -1.246471  -0.5484962   0.025",
  het_group =
    "Sign flips   : n_flip = 6 flip groups  (group order 2^5 = 32, 3 reps)",
  het_est =
    "  log_dist     OLS estimate = -0.8985   95% IPT CI [-1.16, -0.6251]",
  het_p =
    "H0: beta = 0    p-value = 0.031",
  check_alternative =
    "  ? design = \"dyadic_het\" runs the sign-flip test instead: valid under",
  check_alternative_2 =
    "    cells and symmetric about zero, but NOT under additive cluster"
)

## ---- 1. the package still prints what the README claims -------------------
## The exact calls the README shows, verbatim apart from verbose = FALSE (the
## README displays the dispatch banner as a separate block).
fit_qs <- mwperm(y = "log_trade", d = "log_dist",
                 x = c("log_gdp_i", "log_gdp_j"),
                 index = c("importer", "exporter"), data = trade_dyadic,
                 n_reps = 15, seed = 1, verbose = FALSE)
out_qs <- capture.output(print(fit_qs))
stopifnot(readme_lines[["quick_start_est"]] %in% out_qs,
          readme_lines[["quick_start_p"]] %in% out_qs,
          readme_lines[["quick_start_res"]] %in% out_qs)
## the summary() data frame the README shows below the print block; its
## interval columns carry the exact end points at full precision, which is
## where the 0.2.0 bisection numbers survived unnoticed
out_sum <- capture.output(print(summary(fit_qs)))
stopifnot(readme_lines[["quick_start_summary"]] %in% out_sum)

fit_pan <- mwperm(y = "log_trade", d = "fta", x = c("log_gdp_i", "log_gdp_j"),
                  index = c("importer", "exporter", "year"),
                  data = trade_panel, seed = 1, verbose = FALSE)
out_pan <- capture.output(print(fit_pan))
stopifnot(readme_lines[["panel_est"]] %in% out_pan,
          readme_lines[["panel_p"]] %in% out_pan,
          readme_lines[["panel_res"]] %in% out_pan,
          readme_lines[["panel_res_caveat"]] %in% out_pan)

out_chk <- capture.output(print(
  mwperm_check(index = c("importer", "exporter"), data = trade_dyadic)))
stopifnot(readme_lines[["check_resolution"]] %in% out_chk,
          readme_lines[["check_verdict"]] %in% out_chk,
          readme_lines[["check_alternative"]] %in% out_chk,
          readme_lines[["check_alternative_2"]] %in% out_chk)

## the sign-flip transcript under "Extensions"
fit_het <- with(trade_dyadic,
                mwperm_dyadic_het(y = log_trade, d = log_dist,
                                  x = cbind(log_gdp_i, log_gdp_j),
                                  row = importer, col = exporter,
                                  n_flip = 6, n_reps = 3, seed = 1))
out_het <- capture.output(print(fit_het))
stopifnot(readme_lines[["het_group"]] %in% out_het,
          readme_lines[["het_est"]] %in% out_het,
          readme_lines[["het_p"]] %in% out_het)

## the dispatch banner the README prints above the quick-start output
banner <- msgs_of(mwperm(y = "log_trade", d = "log_dist",
                         x = c("log_gdp_i", "log_gdp_j"),
                         index = c("importer", "exporter"),
                         data = trade_dyadic, n_reps = 1, seed = 1,
                         conf_int = FALSE, verbose = TRUE))
stopifnot(any(grepl(paste0("Detected design: dyadic (2 indices, one ",
                           "observation per cell, complete array)"),
                    banner, fixed = TRUE)))

## ---- 2. README.md actually contains those lines ---------------------------
## Skipped where README.md is not shipped alongside the tests (R CMD check).
## Identify it by content, not by path: tests/ has a README.md of its own, so
## a bare "README.md" resolves to the wrong file when run from there.
is_pkg_readme <- function(f)
  file.exists(f) && any(grepl("^# mwperm$", readLines(f, warn = FALSE)))
cands <- c(file.path("..", "README.md"), "README.md")
hit <- Filter(is_pkg_readme, cands)
if (!length(hit)) {
  cat("test-readme.R: package README.md not found, section 2 skipped\n")
} else {
  txt <- readLines(hit[1L], warn = FALSE)
  missing <- names(readme_lines)[
    !vapply(readme_lines, function(l) any(trimws(txt) == trimws(l)),
            logical(1))]
  if (length(missing))
    stop(sprintf(paste0("README.md no longer shows: %s. The package prints ",
                        "these lines, so the README is stale -- update the ",
                        "fenced output blocks in README.md to match, or fix ",
                        "this file if the package's output changed on ",
                        "purpose."),
                 paste(missing, collapse = ", ")), call. = FALSE)
}

passed("test-readme.R")
