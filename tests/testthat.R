## Runner for the testthat suite in tests/testthat/ (the paper-fidelity
## tests). R CMD check runs this file alongside the base-R scripts in tests/;
## see tests/README.md for how the two suites divide the work.
library(testthat)
library(mwperm)

test_check("mwperm")
