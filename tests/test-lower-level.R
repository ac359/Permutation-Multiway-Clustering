## Runner for the files in tests/lower-level-tests/.
##
## R CMD check executes only the .R files at the TOP level of tests/, so the
## lower-level suite would otherwise ship without ever running. This driver
## sources each of those files in turn, in a fresh environment so that fixtures
## with the same name in different files cannot collide, and lets any
## stopifnot() failure propagate -- a failure down there is a failure here.
##
## To run one lower-level file on its own during development:
##   R CMD INSTALL . && Rscript tests/lower-level-tests/test-permset.R
## (install first: those files reach internals through mwperm:::, which resolves
## against the INSTALLED package, never against the working tree.)
library(mwperm)

dir_lower <- if (dir.exists("lower-level-tests")) "lower-level-tests" else
  file.path("tests", "lower-level-tests")
stopifnot(dir.exists(dir_lower))

files <- sort(list.files(dir_lower, pattern = "^test-.*\\.R$",
                         full.names = TRUE))
stopifnot(length(files) > 0L)

cat("test-lower-level.R: running", length(files), "files from",
    sQuote(dir_lower), "\n")
for (f in files) {
  cat("  -> ", basename(f), "\n", sep = "")
  source(f, local = new.env(parent = globalenv()), chdir = FALSE)
}

cat("test-lower-level.R: all assertions passed\n")
