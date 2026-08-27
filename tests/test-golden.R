## Golden-baseline regression gate.
##
## Runs every exported front end on the shipped data at seed = 1, at default and
## at one non-default configuration each, plus find_bicliques() on the
## diagonal-deleted design, and compares the FULL result objects against the
## stored snapshot. Any drift in pvalue / estimate / conf_int / conf_set / K is
## a defect until proven otherwise: the paper's tables and the README's shown
## output are keyed to these seeded numbers.
##
## To re-derive the snapshot after an intentional, documented change:
##   R CMD INSTALL . && Rscript tests/golden/make_baseline.R
## To see which numbers 0.3.0 moved relative to 0.2.0:
##   Rscript tests/golden/make_baseline.R --check --against=baseline-0.2.0.rds
library(mwperm)

## Locate the snapshot whether run from the package root (Rscript tests/...)
## or from tests/ (R CMD check copies the scripts and runs them there).
cands <- c(file.path("tests", "golden", "baseline.rds"),
           file.path("golden", "baseline.rds"))
rds <- cands[file.exists(cands)]
if (!length(rds)) {
  cat("test-golden.R: snapshot not found, skipped\n")
} else {

script <- sub("baseline\\.rds$", "make_baseline.R", rds[1L])
## make_baseline.R writes or checks depending on its arguments; drive the
## check path in this process so a failure is a test failure.
old_args <- commandArgs
env <- new.env(parent = globalenv())
assign("commandArgs", function(trailingOnly = FALSE) "--check", envir = env)
## run from the package root so the script's relative paths resolve
wd <- getwd()
if (basename(wd) == "tests") setwd("..")
script_path <- file.path("tests", "golden", "make_baseline.R")
out <- utils::capture.output(sys.source(script_path, envir = env))
setwd(wd)

bad <- grep("^DIFFERS|^MISSING", out, value = TRUE)
if (length(bad)) {
  cat(paste(out, collapse = "\n"), "\n")
  stop(sprintf(paste0("golden baseline drifted in %d entr%s. Every seeded ",
                      "number in the README, the vignette and the paper is ",
                      "keyed to these; treat this as a defect and find the ",
                      "cause before re-baselining."),
               length(bad), if (length(bad) == 1L) "y" else "ies"),
       call. = FALSE)
}
stopifnot(any(grepl("baseline reproduces", out, fixed = TRUE)))

cat("test-golden.R: all assertions passed\n")
}
