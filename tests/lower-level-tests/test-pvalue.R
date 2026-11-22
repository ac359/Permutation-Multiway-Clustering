## The minorized randomization p-value (GTW 2026 Eq. 10):
##   p(b) = (1 + #{k : min_j a_j(b) <= b_k(b)}) / (K + 1).
## It must live on the grid {1,...,K+1}/(K+1), be conservative on ties, and
## equal 1 at each per-permutation FWL estimate u_j/M_j (the fact the island
## guard uses).
##
## Lower-level test (tests/lower-level-tests/): it exercises the machinery
## beneath the front ends, mostly through package internals, so it must run
## against a FRESHLY INSTALLED copy -- mwperm::: resolves against the
## installed package, never against a source()d working tree.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

set.seed(3)
n <- 7L
K <- 6L
Kp1 <- K + 1L
g <- expand.grid(i = seq_len(n), j = seq_len(n))
N <- nrow(g)
X <- cbind(1, rnorm(N))
D <- matrix(rnorm(N), N, 1)
y <- rnorm(n)[g$i] + rnorm(n)[g$j] + rnorm(N)
op <- mwperm:::.build_obs_perms(cbind(g$i, g$j),
        list(build_perm_set(n, K, seed = 1), build_perm_set(n, K, seed = 2)))
prep <- mwperm:::.ipt_prepare(y, D, X, op, need_perm_D = TRUE)

## ---- 1. p in [1/(K+1), 1] and always a multiple of 1/(K+1) ----------------
for (b in seq(-5, 5, by = 0.05)) {
  p <- mwperm:::.ipt_eval(prep, b)$pvalue
  stopifnot(p >= 1 / Kp1 - 1e-12, p <= 1 + 1e-12)
  stopifnot(abs(p * Kp1 - round(p * Kp1)) < 1e-9)          # on the grid
}

## ---- 2. p(b) is a step function: it changes only at finitely many b -------
## Over a fine sweep it takes only the attainable grid values, and the number of
## distinct values is small (<= K+1).
ps <- vapply(seq(-8, 8, by = 0.01),
             function(b) mwperm:::.ipt_eval(prep, b)$pvalue, numeric(1))
stopifnot(all(ps %in% (seq_len(Kp1) / Kp1)))              # only grid values
stopifnot(length(unique(ps)) <= Kp1)

## ---- 3. conservative ties: b_k == amin counts toward p (>=, not >) --------
## Hand-crafted prep with an exact tie: at b = 0, a = (2,3,4), b = (2,5,6); amin
## = 2
## and b[1] == 2 exactly. Conservative >= must count it, giving p = (1+3)/4 = 1,
## not (1+2)/4 = 0.75.
prep_tie <- list(u = matrix(c(2, 3, 4), 1, 3), v = matrix(c(2, 5, 6), 1, 3),
                 M = array(0, c(1, 1, 3)), W = array(0, c(1, 1, 3)),
                 K = 3L, Kp1 = 4L, d = 1L, has_perm_D = TRUE)
stopifnot(identical(mwperm:::.ipt_eval(prep_tie, 0)$pvalue, 1))

## ---- 4. p == 1 at each per-permutation FWL estimate b = u_j / M_j ---------
## There a_j = 0, so min_j a_j = 0 <= every (nonnegative) b_k: p = 1 exactly.
## This certifies membership and underpins .invert_ci's island guard.
for (j in seq_len(K)) {
  bj <- prep$u[1L, j] / prep$M[1L, 1L, j]
  ev <- mwperm:::.ipt_eval(prep, bj)
  stopifnot(ev$a[j] < 1e-8, identical(ev$pvalue, 1))
}

## ---- 5. smallest attainable p-value is exactly 1/(K+1) --------------------
## Reachable when no permuted statistic meets the minorized observed one.
pmin <- min(vapply(seq(-30, 30, by = 0.1),
                   function(b) mwperm:::.ipt_eval(prep, b)$pvalue, numeric(1)))
stopifnot(identical(pmin, 1 / Kp1))

## ---- 6. a degenerate permutation slice: where Procedure 1 says p = 1 ------
## If `d` lies in the span of [X | X_k] for one permutation k, that slice's
## residualized d is exactly zero in exact arithmetic, so a_k = b_k = 0 and
## Eq. 10 gives p = 1. Numerically the slice is rounding noise, so
## .ipt_prepare(degenerate = "stop") -- the default the engine uses once it
## has ruled out the GLOBAL case -- refuses to compute it, while
## degenerate = "zero" returns the exact-arithmetic answer. Both behaviours
## and the relative 1e-16 threshold on the squared norm are pinned here; the
## behaviour is deliberate (see the roxygen on `degenerate`), not an accident.
## Construction: D is x permuted by group element 2, so for k = 1 the augmented
## design [1, x, x_g] contains D exactly, while D is not in span(X) globally.
x2 <- rnorm(N)
X2 <- cbind(1, x2)
g2 <- op[[2L]]
D2 <- matrix(x2[g2], N, 1)
stopifnot(sum(qr.resid(qr(X2), D2)^2) > 1e-3 * sum(D2^2))  # identified globally
expect_err(mwperm:::.ipt_prepare(y, D2, X2, op, degenerate = "stop",
                                 design = "test"),
           "permutation 1 of the test design")
pz <- mwperm:::.ipt_prepare(y, D2, X2, op, degenerate = "zero")
stopifnot(pz$u[1L, 1L] == 0, pz$v[1L, 1L] == 0, pz$M[1L, 1L, 1L] == 0,
          pz$W[1L, 1L, 1L] == 0,
          any(pz$M[1L, 1L, -1L] > 0))               # the other slices are live
for (b in c(-3, 0, 0.7, 5))
  stopifnot(identical(mwperm:::.ipt_eval(pz, b)$pvalue, 1))
## the threshold is RELATIVE to ||D||^2: rescaling D changes nothing
expect_err(mwperm:::.ipt_prepare(y, 1e6 * D2, X2, op, degenerate = "stop"),
           "No identifying variation")
expect_err(mwperm:::.ipt_prepare(y, 1e-6 * D2, X2, op, degenerate = "stop"),
           "No identifying variation")
## and a slice with real identifying variation never trips it
stopifnot(all(mwperm:::.ipt_prepare(y, D, X, op, degenerate = "stop")$M > 0))

## ---- exact ties with discrete data -------------------------------------------
## Integer y and d with an intercept only: N * (d - mean d) is an integer
## vector, so every a_k and b_k is an exact integer up to the common factor N
## and Eq. (10) can be evaluated in exact integer arithmetic. The package must
## return that p-value. Before the tie tolerance (.tie_tau() in R/core.R) a tie
## that exact arithmetic counts was lost to rounding in about a third of such
## fits, always toward a smaller p-value, and the p-value changed with the row
## order of the data.
local({
  exact_p_int <- function(y, d, ops) {
    w <- length(y) * d - sum(d)          # N * (d - mean d), exact integers
    a <- abs(sum(w * y))                 # the same for every k when X = 1
    b <- vapply(ops[-1L], function(g) abs(sum(w * y[g])), numeric(1))
    (1 + sum(a <= b)) / length(ops)
  }
  set.seed(202)
  n_fit <- 0L
  for (it in 1:120) {
    n <- sample(8:20, 1)
    m <- sample(8:20, 1)
    K <- min(n, m) - 1L
    g <- expand.grid(i = seq_len(n), j = seq_len(m))
    N <- nrow(g)
    d <- rbinom(N, 1, 0.15)
    y <- rbinom(N, 1, 0.5)
    if (length(unique(d)) < 2L) next
    s <- sample(1e5, 1)
    fit <- mwperm_dyadic(y, d, row = g$i, col = g$j, seed = s, conf_int = FALSE)
    pr <- build_perm_set(n, K, seed = mwperm:::.sub_seed(s, 1L))
    pc <- build_perm_set(m, K, seed = mwperm:::.sub_seed(s, 2L))
    key <- g$i + (g$j - 1L) * n
    ops <- lapply(seq_len(K + 1L), function(k)
      match(pr[[k]][g$i] + (pc[[k]][g$j] - 1L) * n, key))
    stopifnot(isTRUE(all.equal(fit$pvalue, exact_p_int(y, d, ops),
                               tolerance = 1e-12)))
    n_fit <- n_fit + 1L
  }
  stopifnot(n_fit >= 100L)

  ## the same data in another row order is the same data: identical p-value
  set.seed(78)
  for (s in 1:40) {
    n <- 15L
    g <- expand.grid(i = seq_len(n), j = seq_len(n))
    N <- nrow(g)
    d <- rbinom(N, 1, 0.3)
    y <- rbinom(N, 1, 0.5)
    o <- sample(N)
    p1 <- mwperm_dyadic(y, d, row = g$i, col = g$j, seed = s,
                        conf_int = FALSE)$pvalue
    p2 <- mwperm_dyadic(y[o], d[o], row = g$i[o], col = g$j[o], seed = s,
                        conf_int = FALSE)$pvalue
    stopifnot(identical(p1, p2))
  }

  ## a prep object without the scale fields keeps the plain comparison, and a
  ## non-finite cross product still gives NA rather than a number
  stopifnot(identical(mwperm:::.tie_tau(list(d = 1L), 0), 0))
  pr0 <- mwperm:::.ipt_prepare(y, matrix(d), matrix(1, N, 1),
                               list(seq_len(N), rev(seq_len(N))))
  stopifnot(mwperm:::.tie_tau(pr0, 0) > 0,
            mwperm:::.tie_tau(pr0, 2) > mwperm:::.tie_tau(pr0, 0),
            length(mwperm:::.tie_tau(pr0, c(-1, 0, 3))) == 3L)
})

passed("test-pvalue.R")
