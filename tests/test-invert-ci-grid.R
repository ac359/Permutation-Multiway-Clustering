## Regression tests for the explicit-`grid` branch of .invert_ci(): the bug fix
## that replaces union (max over reps) aggregation with the MEDIAN across reps
## (Guo, Toulis & Wang 2026, Remark 1), matching the reported p-value and the
## joint-region path, plus the new grid-edge guard, disconnected flag, and
## NA / grid-hygiene handling. See NEWS 0.2.0 "Bug fix (changes reported
## intervals)". Base-R stopifnot style; fast (T8 is opt-in, see below).
library(mwperm)

## ---- helpers ----------------------------------------------------------------
msg_of <- function(expr)
  tryCatch({ expr; NA_character_ },
           error = function(e) conditionMessage(e))
expect_err <- function(expr, pattern) {
  m <- msg_of(expr)
  stopifnot(!is.na(m), grepl(pattern, m, fixed = TRUE))
}

## Build a faithful per-rep prep_list exactly as mwperm_dyadic would (same
## seed scheme: rep r uses seed + r - 1, dims use .sub_seed(., 1:2)), so the
## aggregation rule is exercised in isolation from the front ends.
build_prep_list <- function(y, D, X, coords, n_row, n_col, K, seed, R) {
  lapply(seq_len(R), function(r) {
    s <- seed + r - 1
    Grow <- build_perm_set(n_row, K, seed = mwperm:::.sub_seed(s, 1L))
    Gcol <- build_perm_set(n_col, K, seed = mwperm:::.sub_seed(s, 2L))
    op <- mwperm:::.build_obs_perms(coords, list(Grow, Gcol))
    mwperm:::.ipt_prepare(y, D, X, op, need_perm_D = TRUE)
  })
}
ci_of <- function(pl, grid, alpha, agg)
  mwperm:::.invert_ci(pl, alpha = alpha, centre = 0, scale = 1,
                      y = c(0, 1), D = matrix(c(0, 1)), grid = grid, agg = agg)
width <- function(ci) if (all(is.finite(ci))) diff(ci) else Inf
med_p <- function(pl, b)
  stats::median(vapply(pl, function(pp) mwperm:::.ipt_eval(pp, b)$pvalue,
                       numeric(1)))

## ---- Fixture A: a borderline dyadic design ----------------------------------
## The effect is tuned so H0: beta = 0 sits right on the rejection boundary:
## a majority of reps reject 0 (median p <= alpha) yet a minority accept it,
## the exact regime where the union interval and the reported decision can
## disagree.
set.seed(11)
nA <- 25L
gA <- expand.grid(i = seq_len(nA), j = seq_len(nA))
NA_ <- nrow(gA)
dA <- rnorm(nA)[gA$i] + rnorm(NA_)
yA <- rnorm(nA)[gA$i] + rnorm(nA)[gA$j] + 0.15 * dA + rnorm(NA_)
DA <- matrix(dA, NA_, 1)
XA <- matrix(1, NA_, 1)
coordsA <- cbind(mwperm:::.dense_id(gA$i), mwperm:::.dense_id(gA$j))
KA <- nA - 1L
alpha <- 0.05
grid <- seq(-3, 3, by = 0.01)         # fine grid, 0.01 spacing, contains 0
plA <- build_prep_list(yA, DA, XA, coordsA, nA, nA, KA, seed = 1L, R = 9L)

## ---- T1: deterministic dominance of the aggregation rules -------------------
## median is the tightest de-randomisation; union (= max over reps, the old
## behaviour) and median2 (= min(1, 2*median), valid under arbitrary rep
## dependence) are both no narrower. Guaranteed because acceptance is monotone
## in the pointwise aggregated p-value and the reported CI is the hull of the
## accepted set: p_median <= p_union and p_median <= p_median2 pointwise imply
## the width ordering. NB median2 is NOT bounded above by union (the pointwise
## claim max >= 2*median is false: e.g. reps (.05,.10) give 2*median = .15 >
## max = .10); here 1/(K+1) = 0.04 > alpha/2 so median2 accepts the whole grid.
w_med <- width(ci_of(plA, grid, alpha, "median"))
w_un  <- width(ci_of(plA, grid, alpha, "union"))
w_m2  <- width(ci_of(plA, grid, alpha, "median2"))
stopifnot(w_med <= w_un + 1e-9,        # the fix never widens vs the old union
          w_med <= w_m2 + 1e-9)        # median is the tightest of the three
stopifnot(w_un < w_m2)                 # on this design median2 is the widest

## ---- T2: test / CI coherence (the regression test that matters) -------------
## Reported decision must agree with CI membership: p <= alpha  iff  beta_null
## is outside the interval. True under the median (both invert the SAME median
## p-value); the old union interval breaks it by inverting the maximum.
fitA <- mwperm_dyadic(yA, dA, row = gA$i, col = gA$j,
                      seed = 1, n_reps = 9L, alpha = alpha, grid = grid)
reject <- fitA$pvalue <= fitA$alpha
step <- 0.01
coherent <- function(ci) {
  outside <- (0 < ci[1]) || (0 > ci[2])
  near    <- min(abs(0 - ci[1]), abs(0 - ci[2])) <= step   # one-grid-step slack
  (reject == outside) || near
}
stopifnot(reject)                                  # this fixture rejects 0
stopifnot(coherent(fitA$conf_int))                 # median: decision == CI
## same permutations, union aggregation: 0 lies INSIDE the union interval even
## though the reported (median) decision rejects it -> incoherent by design.
ci_union <- ci_of(plA, grid, alpha, "union")
stopifnot(ci_union[1] <= 0, ci_union[2] >= 0)      # union keeps 0 (contradiction)
stopifnot(!coherent(ci_union))                     # union: decision != CI

## ---- T3: characterise the old pathology (union grows with n_reps) -----------
## Under the nested seed scheme reps accumulate, so the union of per-rep
## acceptance sets can only grow: union width is monotone non-decreasing in
## n_reps (here strictly). The median has no such obligation. Uses the shipped
## trade data, where the growth is strict.
data(trade_dyadic)
yT <- trade_dyadic$log_trade
DT <- matrix(trade_dyadic$log_dist, ncol = 1)
XT <- cbind(1, trade_dyadic$log_gdp_i, trade_dyadic$log_gdp_j)
riT <- mwperm:::.dense_id(trade_dyadic$importer)
ciT <- mwperm:::.dense_id(trade_dyadic$exporter)
coordsT <- cbind(riT, ciT)
nT <- max(riT); KT <- min(nT, max(ciT)) - 1L
gridT <- seq(-3, 1, by = 0.01)
pl1  <- build_prep_list(yT, DT, XT, coordsT, nT, nT, KT, seed = 1L, R = 1L)
pl9  <- build_prep_list(yT, DT, XT, coordsT, nT, nT, KT, seed = 1L, R = 9L)
pl25 <- build_prep_list(yT, DT, XT, coordsT, nT, nT, KT, seed = 1L, R = 25L)
wu1  <- width(ci_of(pl1,  gridT, 0.05, "union"))
wu9  <- width(ci_of(pl9,  gridT, 0.05, "union"))
wu25 <- width(ci_of(pl25, gridT, 0.05, "union"))
stopifnot(wu9 >= wu1 - 1e-9, wu25 >= wu9 - 1e-9)   # monotone non-decreasing
stopifnot(wu9 > wu1)                               # strictly, on this data
## the median interval does not inflate with n_reps
wm1 <- width(ci_of(pl1, gridT, 0.05, "median"))
wm9 <- width(ci_of(pl9, gridT, 0.05, "median"))
stopifnot(wm9 <= wu9)                              # median never wider than union

## ---- T4: grid-median agrees with default bracketing -------------------------
## The two median paths (hull of {b : median_r p_r(b) > alpha} vs median of
## per-rep bracketed end points) coincide when each rep's acceptance set is a
## single interval; assert agreement to within two grid steps.
ci_grid_med <- ci_of(plA, grid, alpha, "median")
if (!isTRUE(attr(ci_grid_med, "disconnected"))) {
  ref <- mwperm:::.ols_reference(yA, DA, XA)
  ci_brack <- mwperm:::.invert_ci(plA, alpha = alpha, centre = ref$estimate,
                                  scale = ref$se, y = yA, D = DA)  # grid = NULL
  stopifnot(!isTRUE(attr(ci_brack, "disconnected")))
  stopifnot(max(abs(as.numeric(ci_grid_med) - as.numeric(ci_brack))) <=
              2 * step + 1e-9)
}

## ---- T5: grid-edge guard (never a silently finite truncated limit) ----------
## A grid too narrow on the upper side: the acceptance region runs off the top
## edge, so the upper limit must be reported as +Inf (the grid does not certify
## a finite bound), the lower limit stays finite, and the truncation is flagged.
narrow <- seq(-0.5, 0.2, by = 0.01)                # acceptance ~ [0.02, 0.46]
ci_n <- ci_of(plA, narrow, alpha, "median")
stopifnot(is.finite(ci_n[1]), ci_n[1] > 0,         # lower end is interior
          is.infinite(ci_n[2]), ci_n[2] > 0)       # upper end clipped -> +Inf
tr <- attr(ci_n, "truncated")
stopifnot(identical(tr, c(FALSE, TRUE)))
gl <- attr(ci_n, "grid_limit")
stopifnot(is.finite(gl[2]), abs(gl[2] - max(narrow)) < 1e-9)   # finite fallback kept
## the engine surfaces this as a note on the fitted object
fit_n <- mwperm_dyadic(yA, dA, row = gA$i, col = gA$j, seed = 1, n_reps = 9L,
                       alpha = alpha, grid = narrow)
stopifnot(any(grepl("reaches the", fit_n$note)),
          any(grepl("grid", fit_n$note)),
          is.infinite(fit_n$conf_int[2]))

## ---- T6: NA safety (a degenerate rep must not corrupt the interval) ---------
K6 <- 24L; Kp1 <- 25L
good6 <- build_prep_list(yA, DA, XA, coordsA, nA, nA, KA, seed = 1L, R = 1L)[[1]]
## all-NaN prep (fully degenerate): p is NA at every b.
bad_all <- list(u = matrix(NaN, 1, K6), v = matrix(0, 1, K6),
                M = array(NaN, c(1, 1, K6)), W = array(0, c(1, 1, K6)),
                K = K6, Kp1 = Kp1, d = 1L, has_perm_D = TRUE)
stopifnot(is.na(mwperm:::.ipt_eval(bad_all, 0)$pvalue))
ci_bad1 <- ci_of(list(bad_all), grid, alpha, "median")
stopifnot(length(ci_bad1) == 2L,
          is.na(ci_bad1[1]) == is.na(ci_bad1[2]))   # never half-NA; no error
## isolated NaN (M = Inf gives Inf*0 = NaN at b = 0 only): mixed with good reps
## the single NA grid point is dropped, the interval stays finite and NA-free.
bad_iso <- list(u = matrix(0, 1, K6), v = matrix(0, 1, K6),
                M = array(Inf, c(1, 1, K6)), W = array(0.5, c(1, 1, K6)),
                K = K6, Kp1 = Kp1, d = 1L, has_perm_D = TRUE)
stopifnot(is.na(mwperm:::.ipt_eval(bad_iso, 0)$pvalue),
          is.finite(mwperm:::.ipt_eval(bad_iso, 0.5)$pvalue))
ci_bad2 <- ci_of(list(good6, bad_iso), grid, alpha, "median")
stopifnot(!anyNA(ci_bad2))                          # no NA endpoint, no error

## ---- T7: grid hygiene (sort / unique / drop non-finite; length-1 errors) ----
clean  <- seq(-1, 1, by = 0.05)
messy  <- c(Inf, sample(c(clean, clean, -Inf, NaN)))   # unsorted + dup + non-finite
stopifnot(identical(as.numeric(ci_of(plA, clean, alpha, "median")),
                    as.numeric(ci_of(plA, messy, alpha, "median"))))
expect_err(ci_of(plA, c(0.5), alpha, "median"),
           "at least two distinct finite values")
expect_err(ci_of(plA, c(0.5, 0.5, Inf), alpha, "median"),
           "at least two distinct finite values")

## ---- T8: slow coverage simulation (opt-in) ----------------------------------
## Set MWPERM_SLOW_TESTS=true to run. 500 null replications at nominal 95%:
## the median-aggregated grid interval covers >= 0.94, and is narrower on
## average than the old union interval.
if (identical(Sys.getenv("MWPERM_SLOW_TESTS"), "true")) {
  set.seed(202)
  nS <- 21L
  gS <- expand.grid(i = seq_len(nS), j = seq_len(nS))
  NS <- nrow(gS)
  coordsS <- cbind(mwperm:::.dense_id(gS$i), mwperm:::.dense_id(gS$j))
  KS <- nS - 1L
  gridS <- seq(-2, 2, by = 0.02)
  R <- 9L; B <- 500L
  cov_med <- 0L; cov_un <- 0L; w_med_s <- 0; w_un_s <- 0
  for (b in seq_len(B)) {
    dS <- rnorm(nS)[gS$i] + rnorm(NS)
    yS <- rnorm(nS)[gS$i] + rnorm(nS)[gS$j] + rnorm(NS)   # true beta = 0
    DS <- matrix(dS, NS, 1); XS <- matrix(1, NS, 1)
    pl <- build_prep_list(yS, DS, XS, coordsS, nS, nS, KS, seed = b, R = R)
    cm <- ci_of(pl, gridS, 0.05, "median")
    cu <- ci_of(pl, gridS, 0.05, "union")
    in_ci <- function(ci) (is.infinite(ci[1]) || ci[1] <= 0) &&
                          (is.infinite(ci[2]) || ci[2] >= 0)
    cov_med <- cov_med + in_ci(cm); cov_un <- cov_un + in_ci(cu)
    w_med_s <- w_med_s + width(cm); w_un_s <- w_un_s + width(cu)
  }
  message(sprintf("T8: coverage median=%.3f union=%.3f ; mean width median=%.3f union=%.3f",
                  cov_med / B, cov_un / B, w_med_s / B, w_un_s / B))
  stopifnot(cov_med / B >= 0.94)              # median interval covers at nominal
  stopifnot(w_med_s < w_un_s)                 # and is tighter than the union
} else {
  message("T8 skipped (set MWPERM_SLOW_TESTS=true to run the coverage sim).")
}

cat("test-invert-ci-grid.R: all assertions passed\n")
