## Phase 2.4 -- the minorized randomization p-value (GTW 2026 Eq. 10):
##   p(b) = (1 + #{k : min_j a_j(b) <= b_k(b)}) / (K + 1).
## It must live on the grid {1,...,K+1}/(K+1), be conservative on ties, and
## equal
## 1 at each per-permutation FWL estimate u_j/M_j (the fact the island guard
## uses).
library(mwperm)

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

## ---- 1. p in [1/(K+1), 1] and always a multiple of 1/(K+1) ------------------
for (b in seq(-5, 5, by = 0.05)) {
  p <- mwperm:::.ipt_eval(prep, b)$pvalue
  stopifnot(p >= 1 / Kp1 - 1e-12, p <= 1 + 1e-12)
  stopifnot(abs(p * Kp1 - round(p * Kp1)) < 1e-9)          # on the grid
}

## ---- 2. p(b) is a step function: it changes only at finitely many b ---------
## Over a fine sweep it takes only the attainable grid values, and the number of
## distinct values is small (<= K+1).
ps <- vapply(seq(-8, 8, by = 0.01),
             function(b) mwperm:::.ipt_eval(prep, b)$pvalue, numeric(1))
stopifnot(all(ps %in% (seq_len(Kp1) / Kp1)))              # only grid values
stopifnot(length(unique(ps)) <= Kp1)

## ---- 3. conservative tie handling: b_k == amin counts toward p (>=, not >)
## ---
## Hand-crafted prep with an exact tie: at b = 0, a = (2,3,4), b = (2,5,6); amin
## = 2
## and b[1] == 2 exactly. Conservative >= must count it, giving p = (1+3)/4 = 1,
## not (1+2)/4 = 0.75.
prep_tie <- list(u = matrix(c(2, 3, 4), 1, 3), v = matrix(c(2, 5, 6), 1, 3),
                 M = array(0, c(1, 1, 3)), W = array(0, c(1, 1, 3)),
                 K = 3L, Kp1 = 4L, d = 1L, has_perm_D = TRUE)
stopifnot(identical(mwperm:::.ipt_eval(prep_tie, 0)$pvalue, 1))

## ---- 4. p == 1 at each per-permutation FWL estimate b = u_j / M_j -----------
## There a_j = 0, so min_j a_j = 0 <= every (nonnegative) b_k: p = 1 exactly.
## This certifies membership and underpins .invert_ci's island guard.
for (j in seq_len(K)) {
  bj <- prep$u[1L, j] / prep$M[1L, 1L, j]
  ev <- mwperm:::.ipt_eval(prep, bj)
  stopifnot(ev$a[j] < 1e-8, identical(ev$pvalue, 1))
}

## ---- 5. smallest attainable p-value is exactly 1/(K+1) ----------------------
## Reachable when no permuted statistic meets the minorized observed one.
pmin <- min(vapply(seq(-30, 30, by = 0.1),
                   function(b) mwperm:::.ipt_eval(prep, b)$pvalue, numeric(1)))
stopifnot(identical(pmin, 1 / Kp1))

cat("test-pvalue.R: all assertions passed\n")
