## Invariance, equivariance and parallel identity: properties that must hold
## for EVERY fit, whatever the design.
##
## Each one is exact and cheap, and each catches a whole class of indexing,
## stacking, partialling or scheduling bug that a Monte-Carlo check would miss.
## The tolerances are not arbitrary. On these fixtures the p-value comes out
## BIT-IDENTICAL under every transformation below, so it is asserted with
## identical(); estimates and interval end points agree only to machine
## epsilon, because rescaling or reordering changes the summation order inside
## the QR. A failure at 1e-10 is a real bug, not drift. (Section 7 is the
## exception: a badly scaled nuisance column goes through a QR of an
## ill-conditioned stack, and its end points are held to 1e-6 of the
## interval's width.)
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

## ---- 7. the units and the origin of a nuisance column change nothing -----
## The projector of Procedure 1 depends on col([X | X_k]) only, so rescaling
## a column of `x`, or shifting it (the intercept is in X), cannot move a
## statistic. The cheap Gram route of .ipt_prepare() squares the condition
## number of the stack and cuts its spectrum relative to the LARGEST column,
## so without a guard a column in large units (or on a large origin) is
## silently not partialled out and the p-value moves. A conditioning guard
## now sends such a stack to a QR of the stack itself (.residualize()). The
## p-value is asserted with identical(); the end points to 1e-6 of the
## interval's width, the precision a badly scaled design can still deliver.
## The fixture has a nuisance column `uw` that is correlated with nothing
## but y, so dropping it from the projection moves the interval (and, for
## some of the transformations below, the p-value).
set.seed(1)
nu <- 20L
gu <- expand.grid(i = seq_len(nu), j = seq_len(nu))
Nu <- nrow(gu)
uz <- runif(nu, 0, 2)[gu$i]
uw <- exp(rnorm(Nu))
ud <- 0.8 * uz + rnorm(Nu)
uy <- 1 + 2 * uz + 0.5 * uw + rnorm(nu)[gu$i] + rnorm(Nu)
## count the calls of the QR route, to show which branch each fit took
ns <- asNamespace("mwperm")
orig_resid <- get(".residualize", envir = ns)
n_qr <- 0L
unlockBinding(".residualize", ns)
assign(".residualize", function(M, V) {
  n_qr <<- n_qr + 1L
  orig_resid(M, V)
}, envir = ns)
lockBinding(".residualize", ns)
qr_calls <- function(expr) {
  n_qr <<- 0L
  force(expr)
  n_qr
}
same_test <- function(f, f_ref)
  identical(f$pvalue, f_ref$pvalue) &&
    identical(f$pvalues_rep, f_ref$pvalues_rep) &&
    max(abs(f$conf_int - f_ref$conf_int)) < 1e-6 * diff(f_ref$conf_int)

## dyadic: units 1e9 and 1e-9 and an origin of 1e6 take the QR route; a
## well-scaled design does not, so its arithmetic is the historical one
nq_ref <- qr_calls(fu0 <- mwperm_dyadic(uy, ud, cbind(uz, uw), gu$i, gu$j,
                                        seed = 1))
stopifnot(nq_ref == 0L)
for (xu in list(cbind(uz, uw * 1e9), cbind(uz, uw * 1e-9),
                cbind(uz, uw + 1e6), cbind(uz * 1e7, uw))) {
  nq <- qr_calls(fu <- mwperm_dyadic(uy, ud, xu, gu$i, gu$j, seed = 1))
  stopifnot(same_test(fu, fu0), nq > 0L)
}
## sign-flip elements take the same guard
hu0 <- mwperm_dyadic_het(uy, ud, cbind(uz, uw), gu$i, gu$j, seed = 1)
nq <- qr_calls(hu <- mwperm_dyadic_het(uy, ud, cbind(uz, uw * 1e9), gu$i,
                                       gu$j, seed = 1))
stopifnot(same_test(hu, hu0), nq > 0L)
## panel with period dummies in X (time_fe = TRUE): every element fixes the
## period, so the dummies repeat in the stack, which is rank-deficient AND
## badly scaled
set.seed(2)
gp <- expand.grid(i = 1:20, j = 1:20, t = 1:3)
pz <- runif(20)[gp$i]
pw <- exp(rnorm(nrow(gp)))
pd <- pz + rnorm(nrow(gp))
py <- 1 + 2 * pz + pw + 0.3 * gp$t + rnorm(20)[gp$i] + rnorm(nrow(gp))
pu0 <- mwperm_panel(py, pd, cbind(pz, pw), gp$i, gp$j, gp$t, seed = 1)
nq <- qr_calls(pu <- mwperm_panel(py, pd, cbind(pz, pw * 1e9), gp$i, gp$j,
                                  gp$t, seed = 1))
stopifnot(pu0$K == 19L, same_test(pu, pu0), nq > 0L)
## two coefficients of interest: the joint region is searched on a grid
## whose step is the naive OLS standard error. solve() refuses the badly
## scaled W'W that standard error used to come from, the step then fell
## back to 1 and the region came back as one or a handful of points, with
## no note. It now comes from the QR that lm.fit() holds whenever solve()
## fails.
ud2 <- 0.5 * uz + rnorm(Nu)
ru0 <- mwperm_dyadic(uy, cbind(ud, ud2), cbind(uz, uw), gu$i, gu$j, seed = 1)
ru <- mwperm_dyadic(uy, cbind(ud, ud2), cbind(uz, uw * 1e9), gu$i, gu$j,
                    seed = 1)
stopifnot(all(is.finite(ru0$se_naive)), all(is.finite(ru$se_naive)),
          max(abs(ru$se_naive / ru0$se_naive - 1)) < 1e-6,
          identical(ru$pvalue, ru0$pvalue),
          nrow(ru0$conf_region) > 10L,
          nrow(ru$conf_region) == nrow(ru0$conf_region),
          max(abs(ru$conf_box - ru0$conf_box)) <
            1e-6 * max(abs(ru0$conf_box)))
unlockBinding(".residualize", ns)
assign(".residualize", orig_resid, envir = ns)
lockBinding(".residualize", ns)

## ---- 8. a serial fit never looks up the core count ------------------------
## parallel::detectCores() shells out on macOS (several ms a call) and used to
## run 2 + n_reps times per fit, also when n_cores = 1, where the count is not
## needed: a default mwperm_irregular() fit (500 repetitions) spent most of
## its wall clock there. Instrument .n_cores_max() to count its calls.
orig_ncm <- get(".n_cores_max", envir = ns)
n_ncm <- 0L
unlockBinding(".n_cores_max", ns)
assign(".n_cores_max", function() {
  n_ncm <<- n_ncm + 1L
  orig_ncm()
}, envir = ns)
lockBinding(".n_cores_max", ns)
f_ser3 <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5,
                        n_reps = 3L)
r_ser <- internal(".plapply")(as.list(1:4), function(i) i + 1L, n_cores = 1L)
n_serial <- n_ncm
## a parallel request still consults it, to clamp n_cores
f_par3 <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 5,
                        n_reps = 3L, n_cores = 2L)
n_parallel <- n_ncm - n_serial
unlockBinding(".n_cores_max", ns)
assign(".n_cores_max", orig_ncm, envir = ns)
lockBinding(".n_cores_max", ns)
stopifnot(n_serial == 0L, n_parallel > 0L,
          identical(r_ser, lapply(as.list(1:4), function(i) i + 1L)),
          isTRUE(same_fit(f_ser3, f_par3)))

## ---- 9. the origin of the covariate of interest changes nothing ----------
## Every design includes the intercept, so the test of d + c is the test of
## d. Far from zero the residualized d is orthogonal to the intercept only
## to rounding, which the origin then multiplies (1e7 sd moved the interval
## by 0.02; from 1e8 sd the fit called d constant and returned p = 1), so
## the engine centres a d whose mean exceeds 1e3 of its standard deviations.
f9 <- mwperm_dyadic(y, d, x = x, row = g$i, col = g$j, seed = 9)
for (k9 in c(1e5, 1e8)) {
  f9k <- mwperm_dyadic(y, d + k9 * sd(d), x = x, row = g$i, col = g$j,
                       seed = 9)
  stopifnot(identical(f9k$pvalue, f9$pvalue),
            max(abs(f9k$conf_int - f9$conf_int)) <
              1e-6 * max(abs(f9$conf_int)))
}
## ... and a d whose only variation beyond x is the rounding of its own
## origin (d = x + 1e9 sd) still has no variation: a warning and p = 1.
w9 <- character(0)
f9u <- withCallingHandlers(
  mwperm_dyadic(y, x[, 1] + 1e9 * sd(x[, 1]), x = x, row = g$i, col = g$j,
                seed = 9, conf_int = FALSE),
  warning = function(w) {
    w9 <<- c(w9, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
stopifnot(f9u$pvalue == 1, any(grepl("no variation", w9, fixed = TRUE)))

passed("test-equivariance.R")
