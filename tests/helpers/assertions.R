## Assertion helpers shared by every test file. See tests/README.md.
##
## Sourced, not run: this file lives in a subdirectory so that R CMD check --
## which runs only the .R files at the top level of tests/ -- does not treat it
## as a test. Each test file loads it with the two-line idiom
##
##   source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
##          else file.path("tests", "helpers", "assertions.R"))
##
## which resolves from the package root (Rscript tests/test-dyadic.R) and from
## tests/ (where R CMD check runs the copied scripts).
##
## Everything here is base R and deliberately tiny: the suite asserts with
## stopifnot(), not with a testing framework, because tests/ ships in the
## tarball and must add no dependency.

## Run `expr`; return its error message, or NA_character_ if it did not error.
msg_of <- function(expr)
  tryCatch({
    expr
    NA_character_
  }, error = function(e) conditionMessage(e))

## Assert that `expr` fails with a message containing `pattern` (fixed, not
## regex: the messages being pinned contain backticks and brackets).
expect_err <- function(expr, pattern) {
  m <- msg_of(expr)
  stopifnot(!is.na(m), grepl(pattern, m, fixed = TRUE))
  invisible(m)
}

## Run `expr` to completion, returning the warning messages it emitted.
warns_of <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  w
}

## Run `expr` to completion, returning the message() output it emitted.
msgs_of <- function(expr) {
  m <- character(0)
  withCallingHandlers(expr, message = function(cnd) {
    m <<- c(m, conditionMessage(cnd))
    invokeRestart("muffleMessage")
  })
  m
}

## Assert that `expr` warns with a message containing `pattern` (fixed).
expect_warn <- function(expr, pattern) {
  w <- warns_of(expr)
  stopifnot(any(grepl(pattern, w, fixed = TRUE)))
  invisible(w)
}

## Compare two "mwperm" fits field by field. Returns TRUE when every field
## agrees under identical(), otherwise the name of the first field that does
## not -- so a failure says WHICH field moved.
same_fit <- function(a, b, skip = "call") {
  for (nm in setdiff(union(names(a), names(b)), skip))
    if (!identical(a[[nm]], b[[nm]])) return(nm)
  TRUE
}

## Reach a package internal both when mwperm is installed (the namespace, which
## is how R CMD check runs these files) and when R/*.R has been source()d into
## the global environment for interactive development.
internal <- function(nm) {
  if (exists(nm, envir = asNamespace("mwperm"), inherits = FALSE))
    get(nm, envir = asNamespace("mwperm"))
  else get(nm, mode = "function")
}

## Print the file's pass line. Called at the end of every test file.
passed <- function(file) cat(file, ": all assertions passed\n", sep = "")
