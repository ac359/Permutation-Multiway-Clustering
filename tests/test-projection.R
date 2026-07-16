## Phase 2.3 — the FWL / partialling-out step (GTW 2026 Eq. 3, Procedure 1 step 1)
## and the affine-in-beta caching that licenses the CI search. Every check
## compares package internals against an INDEPENDENT from-scratch implementation
## (explicit null-space V_k via qr.Q(complete) + explicit projector) that shares
## no statistic/p-value code with R/. See audit/02_correctness.md for the H9
## near-collinear characterisation (finding F2.1) not asserted here.
library(mwperm)

## --- independent references (no shared code with R/'s statistic path) ---------
## orthonormal complement of col(M) from a rank-revealing QR (V' M = 0)
ref_V <- function(M) {
  qM <- qr(M); qr.Q(qM, complete = TRUE)[, seq.int(qM$rank + 1L, nrow(M)), drop = FALSE]
}
## null-space residual via SVD, an entirely separate route from qr.resid
svd_resid <- function(M, V) {
  s <- svd(M); tol <- max(dim(M)) * .Machine$double.eps * max(s$d)
  r <- sum(s$d > tol); U <- s$u[, seq_len(r), drop = FALSE]
  V - U %*% (t(U) %*% V)
}
## naive O(N^3) Procedure 1 p-value at H0: beta = b (explicit V_k, explicit norms)
oracle_pval <- function(y, D, X, obs_perms, b = 0) {
  D <- as.matrix(D); b <- rep(b, length.out = ncol(D))
  ysh <- y - as.numeric(D %*% b); K <- length(obs_perms) - 1L
  a <- numeric(K); bb <- numeric(K)
  for (k in seq_len(K)) {
    g <- obs_perms[[k + 1L]]; V <- ref_V(cbind(X, X[g, , drop = FALSE]))
    P <- V %*% t(V)
    a[k]  <- sqrt(sum(crossprod(D, P %*% ysh)^2))
    bb[k] <- sqrt(sum(crossprod(D, P %*% ysh[g])^2))
  }
  list(pvalue = (1 + sum(bb >= min(a))) / (K + 1L), a = a, b = bb)
}
dyadic_perms <- function(n, K, s) mwperm:::.build_obs_perms(
  as.matrix(expand.grid(i = seq_len(n), j = seq_len(n))),
  list(build_perm_set(n, K, seed = s), build_perm_set(n, K, seed = s + 100L)))

## ---- 1. Eq. (3): V_k' X = 0 and V_k' X_k = 0 (FWL orthogonality) -------------
## The package residualizes via qr.resid(qr([X|X_g]), .); the residuals must be
## orthogonal to BOTH X and the permuted X_g, for y, permuted y, D and permuted D.
set.seed(1)
n <- 8L; K <- 5L
g <- expand.grid(i = seq_len(n), j = seq_len(n)); N <- nrow(g)
X <- cbind(1, rnorm(n)[g$i], rnorm(N))          # well-conditioned nuisance design
D <- matrix(rnorm(N), N, 1); y <- rnorm(n)[g$i] + rnorm(n)[g$j] + rnorm(N)
op <- dyadic_perms(n, K, 20L)
for (k in seq_len(K)) {
  gk <- op[[k + 1L]]; M <- cbind(X, X[gk, ])
  for (v in list(D, y, y[gk], D[gk, , drop = FALSE])) {
    r <- mwperm:::.residualize(M, v)
    stopifnot(max(abs(crossprod(X, r))) < 1e-10)         # orthogonal to X
    stopifnot(max(abs(crossprod(X[gk, ], r))) < 1e-10)   # orthogonal to X_g
  }
}

## ---- 2. rank-deficient [X|X_g] handled exactly (duplicated intercept etc.) ---
## The stacked design ALWAYS duplicates the intercept, and panel time dummies are
## permutation-invariant too. The pivoted qr.resid must give the SAME residuals
## as the SVD null-space projector, to machine precision, on such designs.
gk <- op[[2L]]; M <- cbind(X, X[gk, ])
stopifnot(qr(M)$rank < ncol(M))                          # genuinely rank-deficient
V2 <- matrix(rnorm(N * 2L), N, 2L)
stopifnot(max(abs(mwperm:::.residualize(M, V2) - svd_resid(M, V2))) < 1e-10)

## ---- 3. affine-in-beta: cached .ipt_eval == oracle re-residualising y - D b --
## This is what lets the CI search sweep many b with no new factorisation.
for (d in 1:2) {
  set.seed(30 + d)
  Dm <- matrix(rnorm(N * d), N, d)
  prep <- mwperm:::.ipt_prepare(y, Dm, X, op, need_perm_D = TRUE)
  for (b in list(0, 0.5, -2, rnorm(d))) {
    got <- mwperm:::.ipt_eval(prep, b); ora <- oracle_pval(y, Dm, X, op, b)
    stopifnot(max(abs(got$a - ora$a)) < 1e-9, max(abs(got$b - ora$b)) < 1e-9)
    stopifnot(identical(got$pvalue, ora$pvalue))
  }
}

## ---- 4. oracle equivalence over many random designs (fast subset) -----------
## Full 200-design study with saved output: audit/scripts/02_oracle_equivalence.R.
## Here a fast subset that ships and runs under R CMD check.
set.seed(2026); ncheck <- 0L
for (t in 1:20) {
  nn <- sample(5:8, 1); dd <- sample(1:2, 1); pp <- sample(1:2, 1)
  gg <- expand.grid(i = seq_len(nn), j = seq_len(nn)); NN <- nrow(gg)
  Xt <- cbind(1, matrix(rnorm(NN * pp), NN, pp))
  Dt <- matrix(if (runif(1) < 0.5) rnorm(NN * dd) else exp(0.5 * rnorm(NN * dd)), NN, dd)
  yt <- rnorm(nn)[gg$i] + rnorm(nn)[gg$j] + rnorm(NN)
  Kt <- sample(2:(nn - 1L), 1)
  opt <- mwperm:::.build_obs_perms(cbind(gg$i, gg$j),
           list(build_perm_set(nn, Kt, seed = t), build_perm_set(nn, Kt, seed = t + 7L)))
  prep <- mwperm:::.ipt_prepare(yt, Dt, Xt, opt, need_perm_D = TRUE)
  for (b in list(0, rnorm(dd))) {
    stopifnot(identical(mwperm:::.ipt_eval(prep, b)$pvalue,
                        oracle_pval(yt, Dt, Xt, opt, b)$pvalue))   # EXACT
    ncheck <- ncheck + 1L
  }
}
stopifnot(ncheck == 40L)

cat("test-projection.R: all assertions passed\n")
