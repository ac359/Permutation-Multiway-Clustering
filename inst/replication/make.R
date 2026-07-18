## make.R -- run the whole mwperm replication suite in order.
##
## Runs each numbered script in a fresh R process (the way a reader would run
## them one at a time), records wall-clock time, and writes sessionInfo.txt so
## the environment behind the outputs is on record. Every script is
## cache-backed and resumable: a killed or re-run job reuses completed batches
## under ./cache (redirect with the MWPERM_REPL_CACHE environment variable), so
## re-running make.R after an interruption only computes what is missing.
##
## Usage (from THIS directory, with mwperm installed):
##   Rscript make.R                 # default sim counts (a few minutes total)
##   MC_N=10000 Rscript make.R      # tighter Monte-Carlo error (much longer)
##
## Outputs land in ./out (one .txt table + .rds summary per script). Compare
## ./out against ./expected to confirm a faithful reproduction (small
## Monte-Carlo drift at N below the reference count is expected; the pass/fail
## verdicts printed by each script should agree).

scripts <- c("01_size.R", "02_power.R", "03_ci_coverage.R",
             "04_negative_control.R", "05_permute.R")
rscript <- file.path(R.home("bin"), "Rscript")
dir.create("out", showWarnings = FALSE)

cat("==== mwperm replication suite ====\n")
cat(format(Sys.time()), "\n")
cat("cache:", Sys.getenv("MWPERM_REPL_CACHE", "cache"),
    "| MC_N:", Sys.getenv("MC_N", "(script defaults)"), "\n\n")

timings <- data.frame(script = scripts, seconds = NA_real_, status = "")
for (i in seq_along(scripts)) {
  cat(sprintf("---- [%d/%d] %s ----\n", i, length(scripts), scripts[i]))
  t0 <- Sys.time()
  st <- system2(rscript, scripts[i], stdout = "", stderr = "")
  timings$seconds[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  timings$status[i] <- if (st == 0L) "ok" else sprintf("FAILED (exit %d)", st)
  cat(sprintf("     %s in %.1fs\n\n", timings$status[i], timings$seconds[i]))
}

cat("==== timings ====\n")
print(timings, row.names = FALSE)
cat(sprintf("total: %.1fs\n\n", sum(timings$seconds)))

writeLines(capture.output(sessionInfo()), "out/sessionInfo.txt")
cat("wrote out/sessionInfo.txt\n")
if (any(timings$status != "ok"))
  quit(status = 1L)
