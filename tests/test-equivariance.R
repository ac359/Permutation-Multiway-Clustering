## Invariance, equivariance and parallel identity: properties that must hold
## for EVERY fit, whatever the design.
##
## Each one is exact and cheap, and each catches a whole class of indexing,
## stacking, partialling or scheduling bug that a Monte-Carlo check would miss.
## The tolerances are not arbitrary. On these fixtures the p-value comes out
## BIT-IDENTICAL under every transformation below, so it is asserted with
## identical(); estimates and interval end points agree only to machine
## epsilon, because rescaling or reordering changes the summation order inside
## the QR. A failure at 1e-10 is a real bug, not drift.
##
## Dispatch identity -- mwperm() reproducing the direct front-end call -- is in
## test-main.R, where the dispatcher lives. See tests/README.md.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

## ---- shared dyadic fixture (n = 21 so a 95% CI exists) --------------------
set.seed(1)
n <- 21L
g <- expand.grid(i = seq_len(n), j = seq_len(n))
N <- nrow(g)
x <- cbind(z = rnorm(N))
d <- rnorm(n)[g$i] + rnorm(N)
y <- rnorm(n)[g$i] + rnorm(n)[g$j] + 0.4 * d + rnorm(N)
f0 <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5)

## ---- 1. row order: shuffling the input rows changes nothing ---------------
o <- sample(N)
f1 <- mwperm_dyadic(y[o], d[o], x = x[o, , drop = FALSE],
                    row = g$i[o], col = g$j[o], seed = 5)
stopifnot(identical(f0$pvalue, f1$pvalue),
          identical(f0$pvalues_rep, f1$pvalues_rep),
          max(abs(f0$estimate - f1$estimate)) < 1e-10,
          max(abs(f0$conf_int - f1$conf_int)) < 1e-10)

## ---- 2. nuisance invariance (FWL): y -> y + X c leaves everything alone ---
f2 <- mwperm_dyadic(y + (3 + 1.7 * x[, 1]), d, x = x, row = g$i, col = g$j,
                    seed = 5)
stopifnot(identical(f0$pvalue, f2$pvalue),
          max(abs(f0$estimate - f2$estimate)) < 1e-10,
          max(abs(f0$conf_int - f2$conf_int)) < 1e-10)

## ---- 3. null shift: p(y, beta_null = b) == p(y - D b, beta_null = 0) ------
fa <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5,
                    beta_null = 0.4, conf_int = FALSE)
fb <- mwperm_dyadic(y - 0.4 * d, d, x = x, row = g$i, col = g$j, seed = 5,
                    beta_null = 0, conf_int = FALSE)
stopifnot(identical(fa$pvalue, fb$pvalue))

## ---- 4. scale, sign and affine equivariance -------------------------------
fc <- mwperm_dyadic(y, 2.5 * d, x = x, row = g$i, col = g$j, seed = 5)
stopifnot(identical(f0$pvalue, fc$pvalue),               # p is scale-free
          abs(fc$estimate * 2.5 - f0$estimate) < 1e-10,  # estimate / c
          max(abs(fc$conf_int * 2.5 - f0$conf_int)) < 1e-10)
fs <- mwperm_dyadic(y, -d, x = x, row = g$i, col = g$j, seed = 5)
stopifnot(identical(f0$pvalue, fs$pvalue),
          abs(fs$estimate + f0$estimate) < 1e-10,        # estimate negates
          max(abs(sort(-fs$conf_int) - f0$conf_int)) < 1e-10)
fy <- mwperm_dyadic(2 * y + 7, d, x = x, row = g$i, col = g$j, seed = 5)
stopifnot(identical(f0$pvalue, fy$pvalue),               # p unchanged at b = 0
          abs(fy$estimate / 2 - f0$estimate) < 1e-10,    # estimate scales by a
          max(abs(fy$conf_int / 2 - f0$conf_int)) < 1e-10)

## ---- 5. parallel == serial, bit for bit -----------------------------------
## A README and NEWS promise. Two axes are parallelised: the rep loop (seeded
## per rep) and, at n_reps = 1, the K factorizations (RNG-free). Both must
## reduce order-independently, so the result cannot depend on n_cores.
f_ser <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5,
                       n_reps = 15L)
f_par <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5,
                       n_reps = 15L, n_cores = 2L)
stopifnot(isTRUE(same_fit(f_ser, f_par)))

gm <- expand.grid(i = 1:12, j = 1:12)
gm <- gm[gm$i != gm$j, ]                          # drop the diagonal
set.seed(5)
gm <- gm[sample(nrow(gm), round(0.9 * nrow(gm))), ]
Nm <- nrow(gm)
dm <- rnorm(Nm)
ym <- rnorm(12)[gm$i] + rnorm(12)[gm$j] + rnorm(Nm)
m_ser <- mwperm_missing(ym, dm, row = gm$i, col = gm$j, min_block = 3, seed = 7)
m_par <- mwperm_missing(ym, dm, row = gm$i, col = gm$j, min_block = 3, seed = 7,
                        n_cores = 2L)
stopifnot(isTRUE(same_fit(m_ser, m_par)))

## The three builders that are not plain gather vectors: the sign-flip group
## (signed gathers), the irregular design (a fresh subsample per repetition,
## returned through attr(, "rows")) and the incomplete panel with `L0` (the
## common period set). Rep axis, so n_reps >= 2 with a seed.
h_ser <- mwperm_dyadic_het(y, d, x = x, row = g$i, col = g$j, seed = 5,
                           n_reps = 3L)
h_par <- mwperm_dyadic_het(y, d, x = x, row = g$i, col = g$j, seed = 5,
                           n_reps = 3L, n_cores = 2L)
stopifnot(isTRUE(same_fit(h_ser, h_par)))

set.seed(8)
gi <- expand.grid(i = 1:8, j = 1:8)
gi <- gi[rep(seq_len(nrow(gi)), sample(c(0L, 2L, 3L, 4L), nrow(gi), TRUE)), ]
gi$l <- ave(gi$i, gi$i, gi$j, FUN = seq_along)
Ni <- nrow(gi)
di <- rnorm(64L)[(gi$j - 1L) * 8L + gi$i]           # constant within a cell
yi <- rnorm(8)[gi$i] + rnorm(8)[gi$j] + rnorm(Ni)
i_ser <- mwperm_irregular(yi, di, row = gi$i, col = gi$j, rep = gi$l,
                          L0 = 2L, min_block = 2L, n_reps = 3L, seed = 2)
i_par <- mwperm_irregular(yi, di, row = gi$i, col = gi$j, rep = gi$l,
                          L0 = 2L, min_block = 2L, n_reps = 3L, seed = 2,
                          n_cores = 2L)
stopifnot(isTRUE(same_fit(i_ser, i_par)))

gq <- expand.grid(i = 1:10, j = 1:10, t = 1:3)
set.seed(9)
gq <- gq[!(gq$t == 3L & runif(nrow(gq)) < 0.2), ]  # some pairs miss period 3
Nq <- nrow(gq)
dq <- rnorm(Nq)
yq <- rnorm(10)[gq$i] + rnorm(10)[gq$j] + gq$t + rnorm(Nq)
q_ser <- mwperm_panel_missing(yq, dq, row = gq$i, col = gq$j, time = gq$t,
                              L0 = 2L, n_reps = 3L, seed = 4)
q_par <- mwperm_panel_missing(yq, dq, row = gq$i, col = gq$j, time = gq$t,
                              L0 = 2L, n_reps = 3L, seed = 4, n_cores = 2L)
stopifnot(isTRUE(same_fit(q_ser, q_par)), !is.null(q_ser$periods_used))

## The PSOCK branch (what "auto" selects on Windows) equals serial for a pure
## task, and a worker error propagates on both branches rather than silently
## returning NULL. NOTE: FUN must be an inline anonymous function here --
## .plapply serialises it as an unforced promise, so a bare NAME bound in this
## script's global env is unresolvable on a PSOCK worker. That is a property of
## the test harness, not a package defect: the engine's own calls pass
## functions from the package frame, which every worker can resolve.
X10 <- as.list(1:10)
plapply <- internal(".plapply")
stopifnot(identical(plapply(X10, function(i) i^2 + 1, n_cores = 2L,
                            method = "psock"),
                    lapply(X10, function(i) i^2 + 1)))
## `method = "fork"` is Unix-only: on Windows mclapply() rejects mc.cores > 1
## outright, before the task runs, so the forked branch cannot be exercised
## there. "psock" is exercised everywhere.
err_methods <- if (.Platform$OS.type == "unix") c("fork", "psock") else "psock"
for (m in err_methods) {
  ## suppressWarnings: mclapply emits an expected "core encountered error"
  ## warning before .plapply re-throws the worker error being tested for
  msg <- msg_of(suppressWarnings(
    plapply(as.list(1:4), function(i) if (i == 3L) stop("boom") else i,
            n_cores = 2L, method = m)))
  stopifnot(!is.na(msg), grepl("boom", msg))
}

## a pre-made cluster is reused, not stopped: the engine creates ONE PSOCK
## cluster per fit rather than one per rep on non-fork platforms
cl_pre <- parallel::makePSOCKcluster(2L)
r1 <- plapply(X10, function(i) i^2 + 1, n_cores = 2L, cl = cl_pre)
r2 <- plapply(X10, function(i) i * 2, n_cores = 2L, cl = cl_pre)
stopifnot(identical(r1, lapply(X10, function(i) i^2 + 1)),
          identical(r2, lapply(X10, function(i) i * 2)))   # still usable
parallel::stopCluster(cl_pre)

## n_cores is validated, and clamped to the available cores: silently inside
## .plapply, with a warning naming `n_cores` at the engine entry. mc.cores is
## pinned to 2 so the clamp target is 2 and nothing here ever spawns more than
## R CMD check's two-worker limit.
expect_err(plapply(as.list(1:3), identity, n_cores = NA), "`n_cores`")
old_mc <- getOption("mc.cores")
options(mc.cores = 2L)
stopifnot(identical(plapply(X10, function(i) i + 1L, n_cores = 9999L),
                    lapply(X10, function(i) i + 1L)))
w_clamp <- warns_of(
  f_over <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5,
                          n_reps = 15L, n_cores = 9999L))
options(mc.cores = old_mc)
stopifnot(any(grepl("`n_cores`", w_clamp, fixed = TRUE)),
          isTRUE(same_fit(f_ser, f_over)))          # the clamped run == serial

## ---- 6. seed = NULL forces the rep loop serial ----------------------------
## The RNG-duplication guard. With no seed, forked workers would inherit and
## clone the RNG state, silently drawing the same permutations in every rep --
## a result that looks like n_reps draws but is one. Instrument .plapply in the
## namespace to record how each loop was scheduled: the rep-axis call (length
## n_reps) MUST get n_cores = 1; the RNG-free K-axis calls may use 2.
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
rep_axis <- sched[sched[, "len"] == 3L, , drop = FALSE]      # the rep loop
k_axis <- sched[sched[, "len"] == f0$K, , drop = FALSE]      # the K loops
stopifnot(nrow(rep_axis) >= 1L, all(rep_axis[, "n_cores"] == 1L),
          nrow(k_axis) >= 3L, all(k_axis[, "n_cores"] == 2L),
          length(fit_null$pvalues_rep) == 3L)

passed("test-equivariance.R")
