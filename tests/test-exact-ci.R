## Exact confidence set (Procedure 1, step 3).
##
## The set CI = {b : pval(b) > alpha} is computed in closed form rather than by
## bracketing and bisection: the aggregated p-value is a step function of b
## whose jumps are the roots of |v_k - W_k b| = |u_j - M_j b|, so evaluating it
## at every root and inside every cell between roots determines the set exactly.
## This file pins that against a brute-force evaluation of pval(b) itself.
##
## What is reported is the CLOSURE of {b : pval(b) > alpha}, not the set. An
## acceptance region usually begins and ends strictly inside an open cell of the
## step function, where no boundary value is attained; the exact path then
## reports the bounding root, which the region itself excluded. So a component
## end point may be REJECTED by the test while every interior point is accepted.
## Section 6 pins that relation exactly -- interior accepted, nothing accepted
## outside the hull, every end point either accepted or bounding an accepted
## open cell -- and it is deliberately not the stronger "in the set iff not
## rejected", which is false at the boundary (see NEWS 0.3.0 and CLAUDE.md
## 5.1b). Section 1's evenly spaced brute force cannot see this: its points
## never land on a root.
##
## Reaches internals via ::: -- run against a FRESHLY INSTALLED package.
library(mwperm)

ok <- function(...) stopifnot(...)

## ---- 1. a small dyadic design ---------------------------------------------
set.seed(4)
n <- 8L
g <- expand.grid(i = seq_len(n), j = seq_len(n))
N <- nrow(g)
eta <- rnorm(n); xi <- rnorm(n)
x1 <- rnorm(N)
d <- rnorm(N)
y <- 0.5 * d + 0.3 * x1 + eta[g$i] + xi[g$j] + rnorm(N)

alpha <- 0.25                       # K = 7 => resolution 1/8, so alpha must
                                    # exceed 1/8 for the set to be non-trivial
for (nr in c(1L, 3L)) {
  fit <- mwperm_dyadic(y, d, x = x1, row = g$i, col = g$j,
                       alpha = alpha, n_reps = nr, seed = 7)
  ok(identical(fit$ci_method, "exact"))
  cs <- fit$conf_set
  ok(is.matrix(cs), ncol(cs) == 2L, nrow(cs) >= 1L)
  ## components are ordered and disjoint
  ok(all(cs[, 1L] <= cs[, 2L]))
  if (nrow(cs) > 1L) ok(all(cs[-1L, 1L] > cs[-nrow(cs), 2L]))
  ## conf_int is exactly the hull
  ok(identical(fit$conf_int, c(min(cs[, 1L]), max(cs[, 2L]))))

  ## ---- brute force: pval(b) evaluated directly at each candidate ----------
  ## A coarse check on points that avoid the boundary: away from the roots of
  ## the step function, membership and acceptance coincide. The boundary itself
  ## is section 6's job.
  span <- range(cs[is.finite(cs)])
  pad <- max(diff(span), 1) * 0.25
  ## 201 points: enough to sit inside every cell of the step function on
  ## this design, and each one costs a refit, so the suite stays fast.
  bs <- seq(span[1L] - pad, span[2L] + pad, length.out = 201L)
  inside <- vapply(bs, function(b)
    any(b >= cs[, 1L] & b <= cs[, 2L]), logical(1))
  accepted <- vapply(bs, function(b) {
    pv <- vapply(seq_len(nr), function(r) {
      f <- mwperm_dyadic(y, d, x = x1, row = g$i, col = g$j, beta_null = b,
                         conf_int = FALSE, n_reps = 1L, seed = 7 + r - 1L)
      f$pvalue
    }, numeric(1))
    stats::median(pv) > alpha
  }, logical(1))
  ok(identical(inside, accepted))
}

## ---- 2. the bisection fallback finds the same set --------------------------
## Forcing the fallback (budget 0) must reproduce the exact end points to within
## the bisection tolerance -- the two routes compute the same object.
fit_e <- mwperm_dyadic(y, d, x = x1, row = g$i, col = g$j,
                       alpha = alpha, n_reps = 2L, seed = 7)
old_opt <- options(mwperm.ci_exact_budget = 0)
fit_b <- mwperm_dyadic(y, d, x = x1, row = g$i, col = g$j,
                       alpha = alpha, n_reps = 2L, seed = 7)
options(old_opt)
ok(identical(fit_b$ci_method, "bisection"),
   any(grepl("bracketing and bisection", fit_b$note)),
   isTRUE(all.equal(fit_e$conf_int, fit_b$conf_int, tolerance = 1e-3)))

## ---- 3. the explicit-grid route agrees to the grid spacing -----------------
gr <- seq(fit_e$conf_int[1L] - 0.5, fit_e$conf_int[2L] + 0.5, by = 0.002)
fit_g <- mwperm_dyadic(y, d, x = x1, row = g$i, col = g$j,
                       alpha = alpha, n_reps = 2L, seed = 7, grid = gr)
ok(identical(fit_g$ci_method, "grid"),
   max(abs(fit_g$conf_int - fit_e$conf_int)) <= 0.002 + 1e-12)

## ---- 4. test and interval never disagree away from the boundary -------------
## The whole point of driving both from one rule: a value strictly inside the
## set is not rejected, a value outside it is. (Strictly inside: the end points
## themselves are the closure and are checked in section 6.)
for (b in c(fit_e$conf_int[1L] + 1e-6, mean(fit_e$conf_int),
            fit_e$conf_int[2L] - 1e-6)) {
  f <- mwperm_dyadic(y, d, x = x1, row = g$i, col = g$j, beta_null = b,
                     alpha = alpha, conf_int = FALSE, n_reps = 2L, seed = 7)
  ok(f$pvalue > alpha)
}
for (b in c(fit_e$conf_int[1L] - 0.05, fit_e$conf_int[2L] + 0.05)) {
  f <- mwperm_dyadic(y, d, x = x1, row = g$i, col = g$j, beta_null = b,
                     alpha = alpha, conf_int = FALSE, n_reps = 2L, seed = 7)
  ok(f$pvalue <= alpha)
}

## ---- 5. the breakpoint formulas ---------------------------------------------
## Every root returned must actually solve |v_k - W_k b| = |u_j - M_j b|.
G <- build_perm_set(n, 7L, seed = 11)
op <- mwperm:::.build_obs_perms(cbind(g$i, g$j), list(G, G))
prep <- mwperm:::.ipt_prepare(y, as.matrix(d), cbind(1, x1), op)
roots <- mwperm:::.ci_breakpoints(list(prep))
u <- prep$u[1L, ]; M <- prep$M[1L, 1L, ]
v <- prep$v[1L, ]; W <- prep$W[1L, 1L, ]
ok(length(roots) > 0L, all(is.finite(roots)), !is.unsorted(roots),
   !anyDuplicated(roots))
hits <- vapply(roots, function(b) {
  a <- abs(u - M * b); bk <- abs(v - W * b)
  min(abs(outer(a, bk, `-`))) <= 1e-8 * max(1, max(a), max(bk))
}, logical(1))
ok(all(hits))

## ---- 6. the boundary: what is reported is the CLOSURE of the set ------------
## Sections 1 and 4 probe points that are never roots of the step function, so
## neither can see the boundary. Here the roots themselves are probed, on the
## fit's OWN permutations (rebuilt from the documented seed scheme: rep r uses
## seed + r - 1, and dimension j uses .sub_seed(rep_seed, j)). The relation
## asserted is the closure relation, NOT "in the set iff not rejected" -- the
## latter is false at an end point that bounds an accepted open cell.
mk_preps <- function(nr, K) lapply(seq_len(nr), function(r) {
  rs <- 7L + r - 1L                  # the fit was run with seed = 7
  Grow <- build_perm_set(n, K, seed = mwperm:::.sub_seed(rs, 1L))
  Gcol <- build_perm_set(n, K, seed = mwperm:::.sub_seed(rs, 2L))
  mwperm:::.ipt_prepare(y, as.matrix(d), cbind(1, x1),
                        mwperm:::.build_obs_perms(cbind(g$i, g$j),
                                                  list(Grow, Gcol)),
                        need_perm_D = TRUE)
})

n_rejected_end <- 0L                 # end points the test itself rejects
for (nr in c(1L, 3L)) {
  fit <- mwperm_dyadic(y, d, x = x1, row = g$i, col = g$j,
                       alpha = alpha, n_reps = nr, seed = 7)
  pl <- mk_preps(nr, fit$K)
  ## p-value by the per-rep route (.ipt_eval), independent of the vectorized
  ## .pval_matrix the exact path uses.
  pval_at <- function(b)
    stats::median(vapply(pl, function(pp) mwperm:::.ipt_eval(pp, b)$pvalue,
                         numeric(1)))
  ## the rebuilt permutations ARE the fit's: same p-value at the null
  ok(identical(pval_at(fit$beta_null), fit$pvalue))

  roots <- mwperm:::.ci_breakpoints(pl)
  m <- length(roots)
  ok(m >= 2L)
  ## one atom per root, one interior point per cell between roots, one beyond
  ## each end: the p-value is constant on each cell, so these exhaust the line.
  mids  <- (roots[-m] + roots[-1L]) / 2
  atoms <- sort(c(roots[1L] - (1 + abs(roots[1L])), roots, mids,
                  roots[m] + (1 + abs(roots[m]))))
  acc <- vapply(atoms, function(b) pval_at(b) > alpha, logical(1))

  cs <- fit$conf_set
  hull <- fit$conf_int
  ## (a) conservative: nothing the test accepts lies outside the reported hull
  ok(!any(acc & (atoms < hull[1L] | atoms > hull[2L])))

  for (i in seq_len(nrow(cs))) {
    lo <- cs[i, 1L]; hi <- cs[i, 2L]
    ## (b) every point strictly inside a component is accepted
    ok(all(acc[atoms > lo & atoms < hi]))
    ## (c) each finite end point is a root, and is either accepted itself or is
    ## the infimum / supremum of an accepted open cell (the closure case)
    if (is.finite(lo)) {
      ok(any(roots == lo))
      j <- match(lo, atoms)
      in_lo <- pval_at(lo) > alpha
      ok(in_lo || (j < length(atoms) && acc[j + 1L]))
      if (!in_lo) n_rejected_end <- n_rejected_end + 1L
    }
    if (is.finite(hi)) {
      ok(any(roots == hi))
      j <- match(hi, atoms)
      in_hi <- pval_at(hi) > alpha
      ok(in_hi || (j > 1L && acc[j - 1L]))
      if (!in_hi) n_rejected_end <- n_rejected_end + 1L
    }
  }
}
## The closure case is REACHED on this fixture, so section 6 is not vacuous and
## the stronger iff claim is demonstrably false here: at least one reported end
## point is a value the test rejects. If this ever fails because the exact path
## moved to the attained side, that is a numeric change -- update NEWS, the
## golden baseline and the docs rather than deleting the assertion.
ok(n_rejected_end > 0L)

cat("test-exact-ci.R: all assertions passed\n")
