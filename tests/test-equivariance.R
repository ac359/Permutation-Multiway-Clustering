## Equivariance / invariance / identity contracts. All exact, all
## cheap, all high-yield: each catches a whole class of indexing, stacking,
## partialling or scheduling bug that Monte Carlo would miss.
##
## Tolerances are not arbitrary. On these fixtures the p-value comes out
## BIT-IDENTICAL under every transformation below, so it is asserted with
## identical(). Estimates and interval endpoints agree only to machine
## epsilon, because rescaling or reordering changes the summation order
## inside the QR; those are asserted at 1e-10. A failure at 1e-10 is a real
## bug, not drift.
library(mwperm)

same_fit <- function(a, b, skip = "call") {
  for (nm in setdiff(union(names(a), names(b)), skip))
    if (!identical(a[[nm]], b[[nm]])) return(nm)
  TRUE
}

## ---- shared dyadic fixture (n = 21 so a 95% CI exists) ----------------------
set.seed(1)
n <- 21L
g <- expand.grid(i = seq_len(n), j = seq_len(n))
N <- nrow(g)
x <- cbind(z = rnorm(N))
d <- rnorm(n)[g$i] + rnorm(N)
y <- rnorm(n)[g$i] + rnorm(n)[g$j] + 0.4 * d + rnorm(N)
f0 <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5)

## ---- 1. row-order invariance: shuffling input rows changes nothing ----------
o <- sample(N)
f1 <- mwperm_dyadic(y[o], d[o], x = x[o, , drop = FALSE],
                    row = g$i[o], col = g$j[o], seed = 5)
stopifnot(identical(f0$pvalue, f1$pvalue))               # bit-identical p
stopifnot(identical(f0$pvalues_rep, f1$pvalues_rep))
stopifnot(max(abs(f0$estimate - f1$estimate)) < 1e-10)   # OLS to machine eps
stopifnot(max(abs(f0$conf_int - f1$conf_int)) < 1e-10)

## ---- 2. nuisance invariance (FWL): y -> y + X c leaves everything unchanged
## --
f2 <- mwperm_dyadic(y + (3 + 1.7 * x[, 1]), d, x = x, row = g$i, col = g$j,
                    seed = 5)
stopifnot(identical(f0$pvalue, f2$pvalue))
stopifnot(max(abs(f0$estimate - f2$estimate)) < 1e-10)
stopifnot(max(abs(f0$conf_int - f2$conf_int)) < 1e-10)

## ---- 3. null shift: p(y, beta_null = b) == p(y - D b, beta_null = 0) --------
fa <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5,
                    beta_null = 0.4, conf_int = FALSE)
fb <- mwperm_dyadic(y - 0.4 * d, d, x = x, row = g$i, col = g$j, seed = 5,
                    beta_null = 0, conf_int = FALSE)
stopifnot(identical(fa$pvalue, fb$pvalue))

## ---- 4. scale / sign / affine equivariance ----------------------------------
fc <- mwperm_dyadic(y, 2.5 * d, x = x, row = g$i, col = g$j, seed = 5)
stopifnot(identical(f0$pvalue, fc$pvalue))               # p scale-free
stopifnot(abs(fc$estimate * 2.5 - f0$estimate) < 1e-10)  # estimate / c
stopifnot(max(abs(fc$conf_int * 2.5 - f0$conf_int)) < 1e-10)
fs <- mwperm_dyadic(y, -d, x = x, row = g$i, col = g$j, seed = 5)
stopifnot(identical(f0$pvalue, fs$pvalue))
stopifnot(abs(fs$estimate + f0$estimate) < 1e-10)        # estimate negates
stopifnot(max(abs(sort(-fs$conf_int) - f0$conf_int)) < 1e-10)  # CI negates
fy <- mwperm_dyadic(2 * y + 7, d, x = x, row = g$i, col = g$j, seed = 5)
stopifnot(identical(f0$pvalue, fy$pvalue))               # p unchanged at b = 0
stopifnot(abs(fy$estimate / 2 - f0$estimate) < 1e-10)    # estimate scales by a
stopifnot(max(abs(fy$conf_int / 2 - f0$conf_int)) < 1e-10)

## ---- 5. dispatch identity: mwperm() == the direct front-end, all 5 designs
## ---
## (vignette/man/NEWS promise; fields identical except `call`, plus the `auto`
## metadata mwperm() adds)
fd <- mwperm(y = y, d = d, x = x, index = list(row = g$i, col = g$j),
             design = "dyadic", seed = 5, verbose = FALSE)
stopifnot(isTRUE(same_fit(f0, fd, skip = c("call", "auto"))))

gp <- expand.grid(i = 1:6, j = 1:6, t = 1:4)
Np <- nrow(gp)
set.seed(2)
dp <- rnorm(Np)
yp <- rnorm(6)[gp$i] + rnorm(6)[gp$j] + cumsum(rnorm(4))[gp$t] + 0.3 * dp +
  rnorm(Np)
p_dir <- mwperm_panel(yp, dp, row = gp$i, col = gp$j, time = gp$t, seed = 3)
p_dis <- mwperm(y = yp, d = dp, index = list(row = gp$i, col = gp$j),
                time = gp$t, design = "panel", seed = 3, verbose = FALSE)
stopifnot(isTRUE(same_fit(p_dir, p_dis, skip = c("call", "auto"))))

g3 <- expand.grid(a = 1:5, b = 1:5, c = 1:5)
N3 <- nrow(g3)
set.seed(3)
d3 <- rnorm(N3)
y3 <- rnorm(5)[g3$a] + rnorm(5)[g3$b] + rnorm(5)[g3$c] + rnorm(N3)
t_dir <- mwperm_threeway(y3, d3, id1 = g3$a, id2 = g3$b, id3 = g3$c, seed = 4)
t_dis <- mwperm(y = y3, d = d3, index = list(id1 = g3$a, id2 = g3$b,
                                             id3 = g3$c),
                design = "threeway", seed = 4, verbose = FALSE)
stopifnot(isTRUE(same_fit(t_dir, t_dis, skip = c("call", "auto"))))

gl <- expand.grid(l = 1:5, i = 1:4, j = 1:4)
Nl <- nrow(gl)
set.seed(4)
dl <- rnorm(Nl)
yl <- rnorm(16)[as.integer(interaction(gl$i, gl$j))] + 0.5 * dl + rnorm(Nl)
l_dir <- mwperm_layout(yl, dl, row = gl$i, col = gl$j, seed = 6)
l_dis <- suppressWarnings(mwperm(y = yl, d = dl, index = list(row = gl$i,
                                                              col = gl$j),
                                 design = "layout", seed = 6, verbose = FALSE))
stopifnot(isTRUE(same_fit(l_dir, l_dis, skip = c("call", "auto"))))

gm <- expand.grid(i = 1:12, j = 1:12)
gm <- gm[gm$i != gm$j, ]     # drop diagonal
set.seed(5)
keep <- sample(nrow(gm), round(0.9 * nrow(gm)))
gm <- gm[keep, ]
Nm <- nrow(gm)
dm <- rnorm(Nm)
ym <- rnorm(12)[gm$i] + rnorm(12)[gm$j] + rnorm(Nm)
m_dir <- mwperm_missing(ym, dm, row = gm$i, col = gm$j, min_block = 3, seed = 7)
m_dis <- mwperm(y = ym, d = dm, index = list(row = gm$i, col = gm$j),
                design = "missing", min_block = 3, seed = 7, verbose = FALSE)
stopifnot(isTRUE(same_fit(m_dir, m_dis, skip = c("call", "auto"))))

## ---- 6. parallel identity: n_cores = 2 (fork) bit-identical to serial -------
## (README/NEWS promise). Rep axis: seeded n_reps = 15; K axis: n_reps = 1.
f_ser <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5,
                       n_reps = 15L)
f_par <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5,
                       n_reps = 15L,
                       n_cores = 2L)
stopifnot(isTRUE(same_fit(f_ser, f_par)))
m_par <- mwperm_missing(ym, dm, row = gm$i, col = gm$j, min_block = 3, seed = 7,
                        n_cores = 2L)
stopifnot(isTRUE(same_fit(m_dir, m_par)))

## .plapply PSOCK branch (the Windows fallback) equals serial for a pure task,
## and worker errors propagate on both branches. NOTE: FUN must be passed as an
## inline anonymous function here -- .plapply serialises FUN as an unforced
## promise, so a bare *name* bound in the test's global env is unresolvable on a
## PSOCK worker. This is a property of the test harness, not a package defect:
## the engine's own calls pass functions from the package frame, which every
## worker can resolve.
X10 <- as.list(1:10)
stopifnot(identical(
  mwperm:::.plapply(X10, function(i) i^2 + 1, n_cores = 2L, method = "psock"),
  lapply(X10, function(i) i^2 + 1)))
## `method = "fork"` is a Unix-only code path: on Windows parallel::mclapply()
## rejects mc.cores > 1 outright ("'mc.cores' > 1 is not supported on Windows")
## before the task ever runs, so the forked branch cannot be exercised there.
## Test it only where forking exists; "psock" -- the branch .plapply's "auto"
## selects on Windows -- is exercised on every platform.
err_methods <- if (.Platform$OS.type == "unix") c("fork", "psock") else "psock"
for (m in err_methods) {
  ## suppressWarnings: mclapply emits an expected "core encountered error"
  ## warning before .plapply re-throws the worker error we are testing for
  msg <- tryCatch({
    suppressWarnings(mwperm:::.plapply(
      as.list(1:4),
      function(i) if (i == 3L) stop("boom") else i,
      n_cores = 2L, method = m))
    NA_character_
  }, error = function(e) conditionMessage(e))
  stopifnot(!is.na(msg), grepl("boom", msg))
}

## a pre-made cluster is reused, not stopped (F6.1 fix: the engine creates ONE
## PSOCK cluster per fit instead of one per rep on non-fork platforms)
cl_pre <- parallel::makePSOCKcluster(2L)
r1 <- mwperm:::.plapply(X10, function(i) i^2 + 1, n_cores = 2L, cl = cl_pre)
r2 <- mwperm:::.plapply(X10, function(i) i * 2, n_cores = 2L, cl = cl_pre)
stopifnot(identical(r1, lapply(X10, function(i) i^2 + 1)),
          identical(r2, lapply(X10, function(i) i * 2)))  # still usable
parallel::stopCluster(cl_pre)

## n_cores is clamped to the available core count (F6.5 fix): silently inside
## .plapply, with a warning naming `n_cores` at the engine entry. mc.cores is
## pinned to 2 so the clamp target is 2 and the assertions never spawn more
## than R CMD check's 2-worker limit (the clamp also honours mc.cores).
old_mc <- getOption("mc.cores")
options(mc.cores = 2L)
stopifnot(identical(
  mwperm:::.plapply(X10, function(i) i + 1L, n_cores = 9999L),
  lapply(X10, function(i) i + 1L)))
w_clamp <- character(0)
f_over <- withCallingHandlers(
  mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5, n_reps = 15L,
                n_cores = 9999L),
  warning = function(cnd) {
    w_clamp <<- c(w_clamp, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
options(mc.cores = old_mc)
stopifnot(any(grepl("`n_cores`", w_clamp, fixed = TRUE)),
          isTRUE(same_fit(f_ser,
                          f_over)))                 # clamped run == serial

## ---- 7. seed = NULL forces the rep loop serial (the RNG-duplication guard)
## ---
## Instrument .plapply in the namespace to record how each loop is scheduled:
## with seed = NULL and n_cores = 2, the rep-axis call (length n_reps) MUST get
## n_cores = 1 (else forked workers would clone the RNG state and silently
## duplicate permutation draws); the RNG-free K-axis calls may use 2.
calls <- list()
ns <- asNamespace("mwperm")
orig_plapply <- get(".plapply", envir = ns)
wrap <- function(X, FUN, n_cores = 1L, method = c("auto", "fork", "psock"),
                 cl = NULL) {
  calls[[length(calls) + 1L]] <<- c(len = length(X),
                                    n_cores = as.integer(n_cores))
  orig_plapply(X, FUN, n_cores = n_cores, method = method, cl = cl)
}
unlockBinding(".plapply", ns)
assign(".plapply", wrap, envir = ns)
lockBinding(".plapply", ns)
fit_null <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = NULL,
                          n_reps = 3L, n_cores = 2L, conf_int = FALSE)
unlockBinding(".plapply", ns)
assign(".plapply", orig_plapply, envir = ns)
lockBinding(".plapply", ns)
sched <- do.call(rbind, calls)
rep_axis_calls <- sched[sched[, "len"] == 3L, , drop = FALSE]   # the rep loop
stopifnot(nrow(rep_axis_calls) >= 1L,
          all(rep_axis_calls[, "n_cores"] == 1L))               # guard fired
k_axis_calls <- sched[sched[, "len"] == f0$K, , drop = FALSE]   # the K loops
stopifnot(nrow(k_axis_calls) >= 3L,
          all(k_axis_calls[,
                           "n_cores"] == 2L))                 # K axis parallel
stopifnot(length(fit_null$pvalues_rep) == 3L)

cat("test-equivariance.R: all assertions passed\n")
