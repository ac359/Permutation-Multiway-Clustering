## The FWL / partialling-out step (GTW 2026 Eq. 3, Procedure 1 step 1) and
## the affine-in-beta caching that licenses the CI search. Every check
## compares package internals against an INDEPENDENT from-scratch implementation
## (explicit null-space V_k via qr.Q(complete) + explicit projector) that shares
## no statistic/p-value code with R/. Not asserted here: behaviour when `d` is
## near-collinear with `x` but not exactly so. The statistic stays correct, but
## the residualized quantities lose relative precision as the angle closes, so
## any tolerance-based assertion would be a test of the fixture rather than of
## the code. The exactly-degenerate case IS covered, in test-edgecases.R.
##
## Lower-level test (tests/lower-level-tests/): it exercises the machinery
## beneath the front ends, mostly through package internals, so it must run
## against a FRESHLY INSTALLED copy -- mwperm::: resolves against the
## installed package, never against a source()d working tree.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

## --- independent references (no shared code with R/'s statistic path) --------
## orthonormal complement of col(M) from a rank-revealing QR (V' M = 0)
ref_V <- function(M) {
  qM <- qr(M)
  qr.Q(qM, complete = TRUE)[, seq.int(qM$rank + 1L, nrow(M)), drop = FALSE]
}
## null-space residual via SVD, an entirely separate route from qr.resid
svd_resid <- function(M, V) {
  s <- svd(M)
  tol <- max(dim(M)) * .Machine$double.eps * max(s$d)
  r <- sum(s$d > tol)
  U <- s$u[, seq_len(r), drop = FALSE]
  V - U %*% (t(U) %*% V)
}
## naive O(N^3) Procedure 1 p-value at H0: beta = b (explicit V_k, explicit
## norms)
oracle_pval <- function(y, D, X, obs_perms, b = 0) {
  D <- as.matrix(D)
  b <- rep(b, length.out = ncol(D))
  ysh <- y - as.numeric(D %*% b)
  K <- length(obs_perms) - 1L
  a <- numeric(K)
  bb <- numeric(K)
  for (k in seq_len(K)) {
    g <- obs_perms[[k + 1L]]
    V <- ref_V(cbind(X, X[g, , drop = FALSE]))
    P <- V %*% t(V)
    a[k]  <- sqrt(sum(crossprod(D, P %*% ysh)^2))
    bb[k] <- sqrt(sum(crossprod(D, P %*% ysh[g])^2))
  }
  list(pvalue = (1 + sum(bb >= min(a))) / (K + 1L), a = a, b = bb)
}
dyadic_perms <- function(n, K, s) mwperm:::.build_obs_perms(
  as.matrix(expand.grid(i = seq_len(n), j = seq_len(n))),
  list(build_perm_set(n, K, seed = s), build_perm_set(n, K, seed = s + 100L)))

## ---- 1. Eq. (3): V_k' X = 0 and V_k' X_k = 0 (FWL orthogonality) ----------
## .residualize() (qr.resid of the stack [X|X_g]: the fallback of the fit path,
## see section 5) must return residuals orthogonal to BOTH X and the permuted
## X_g, for y, permuted y, D and permuted D.
set.seed(1)
n <- 8L
K <- 5L
g <- expand.grid(i = seq_len(n), j = seq_len(n))
N <- nrow(g)
X <- cbind(1, rnorm(n)[g$i],
           rnorm(N))          # well-conditioned nuisance design
D <- matrix(rnorm(N), N, 1)
y <- rnorm(n)[g$i] + rnorm(n)[g$j] + rnorm(N)
op <- dyadic_perms(n, K, 20L)
for (k in seq_len(K)) {
  gk <- op[[k + 1L]]
  M <- cbind(X, X[gk, ])
  for (v in list(D, y, y[gk], D[gk, , drop = FALSE])) {
    r <- mwperm:::.residualize(M, v)
    stopifnot(max(abs(crossprod(X, r))) < 1e-10)         # orthogonal to X
    stopifnot(max(abs(crossprod(X[gk, ], r))) < 1e-10)   # orthogonal to X_g
  }
}

## ---- 2. rank-deficient [X|X_g] handled exactly ----------------------------
## The stacked design ALWAYS duplicates the intercept, and panel time dummies
## are
## permutation-invariant too. The pivoted qr.resid must give the SAME residuals
## as the SVD null-space projector, to machine precision, on such designs.
gk <- op[[2L]]
M <- cbind(X, X[gk, ])
stopifnot(qr(M)$rank < ncol(M))                # genuinely rank-deficient
V2 <- matrix(rnorm(N * 2L), N, 2L)
stopifnot(max(abs(mwperm:::.residualize(M, V2) - svd_resid(M, V2))) < 1e-10)

## ---- 3. affine-in-beta: cached .ipt_eval == oracle re-residualising y - D b
## --
## This is what lets the CI search sweep many b with no new factorisation.
for (d in 1:2) {
  set.seed(30 + d)
  Dm <- matrix(rnorm(N * d), N, d)
  prep <- mwperm:::.ipt_prepare(y, Dm, X, op, need_perm_D = TRUE)
  for (b in list(0, 0.5, -2, rnorm(d))) {
    got <- mwperm:::.ipt_eval(prep, b)
    ora <- oracle_pval(y, Dm, X, op, b)
    stopifnot(max(abs(got$a - ora$a)) < 1e-9, max(abs(got$b - ora$b)) < 1e-9)
    stopifnot(identical(got$pvalue, ora$pvalue))
  }
}

## ---- 4. oracle equivalence over many random designs (fast subset) ---------
## A fast subset sized to run under R CMD check. The same comparison has been
## run over 200 random designs off-line; nothing here depends on that, and
## widening `ncheck` below reproduces it at proportionate cost.
set.seed(2026)
ncheck <- 0L
for (t in 1:20) {
  nn <- sample(5:8, 1)
  dd <- sample(1:2, 1)
  pp <- sample(1:2, 1)
  gg <- expand.grid(i = seq_len(nn), j = seq_len(nn))
  NN <- nrow(gg)
  Xt <- cbind(1, matrix(rnorm(NN * pp), NN, pp))
  Dt <- matrix(if (runif(1) < 0.5) rnorm(NN * dd)
               else exp(0.5 * rnorm(NN * dd)), NN, dd)
  yt <- rnorm(nn)[gg$i] + rnorm(nn)[gg$j] + rnorm(NN)
  Kt <- sample(2:(nn - 1L), 1)
  opt <- mwperm:::.build_obs_perms(cbind(gg$i, gg$j),
           list(build_perm_set(nn, Kt, seed = t),
                build_perm_set(nn, Kt, seed = t + 7L)))
  prep <- mwperm:::.ipt_prepare(yt, Dt, Xt, opt, need_perm_D = TRUE)
  for (b in list(0, rnorm(dd))) {
    stopifnot(identical(mwperm:::.ipt_eval(prep, b)$pvalue,
                        oracle_pval(yt, Dt, Xt, opt, b)$pvalue))   # EXACT
    ncheck <- ncheck + 1L
  }
}
stopifnot(ncheck == 40L)

## ---- 5. the two routes of .ipt_prepare(): Gram pseudo-inverse and QR ------
## .ipt_prepare() normally residualizes D through an eigendecomposition of the
## 2p x 2p Gram matrix of [X | X_k]; when that matrix is too badly conditioned
## (a nuisance column in large units, or with a large mean and a small spread)
## a guard sends the element to .residualize(), a QR of the stack. Pinned here:
##  (a) on a well-conditioned design the result agrees with the QR route to
##      rounding (that the Gram route is the one taken there is pinned by
##      test-equivariance.R section 7);
##  (b) on the same design with one column in units of 1e9 every element
##      takes the QR route -- its cached cross products are EXACTLY those of
##      .residualize() -- and they agree with (a)'s, so rescaling a column
##      moved nothing. Permutations and sign flips alike.
set.seed(5)
n5 <- 10L
g5 <- expand.grid(i = seq_len(n5), j = seq_len(n5))
N5 <- nrow(g5)
z5 <- rnorm(n5)[g5$i]
w5 <- exp(rnorm(N5))
X5 <- cbind(1, z5, w5)
D5 <- cbind(0.5 * z5 + rnorm(N5), rnorm(N5))         # d = 2
y5 <- 1 + z5 + w5 + rnorm(N5)
fl <- build_flip_set(n5, n5, 4L, seed = 3)
ops5 <- c(dyadic_perms(n5, 4L, 9L),
          lapply(fl[-1L], function(f) list(g = NULL, s = f$row[g5$i] *
                                                         f$col[g5$j])))
## by hand, through the QR route, for every non-identity element
by_qr <- function(X) {
  sl <- lapply(ops5[-1L], function(op) {
    Xg <- mwperm:::.apply_op(op, X)
    Dr <- mwperm:::.residualize(cbind(X, Xg), D5)
    list(u = crossprod(Dr, y5), v = crossprod(Dr, mwperm:::.apply_op(op, y5)),
         M = crossprod(Dr), W = crossprod(Dr, mwperm:::.apply_op(op, D5)))
  })
  list(u = sapply(sl, `[[`, "u"), v = sapply(sl, `[[`, "v"),
       M = simplify2array(lapply(sl, `[[`, "M")),
       W = simplify2array(lapply(sl, `[[`, "W")))
}
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
X5s <- X5
X5s[, 3L] <- X5s[, 3L] * 1e9
pg <- mwperm:::.ipt_prepare(y5, D5, X5, ops5, need_perm_D = TRUE)
ps <- mwperm:::.ipt_prepare(y5, D5, X5s, ops5, need_perm_D = TRUE)
q0 <- by_qr(X5)
qs <- by_qr(X5s)
for (f in c("u", "v", "M", "W")) {
  stopifnot(rel(pg[[f]], q0[[f]]) < 1e-10,                 # (a)
            identical(unname(ps[[f]]), unname(qs[[f]])),     # (b) QR route
            rel(ps[[f]], pg[[f]]) < 1e-8)                    # (b) invariance
}
## the p-value at several nulls is the same object on both designs
for (b in list(c(0, 0), c(0.5, -0.2), c(2, 1)))
  stopifnot(identical(mwperm:::.ipt_eval(ps, b)$pvalue,
                      mwperm:::.ipt_eval(pg, b)$pvalue))

## ---- 6. an element that fixes d and x keeps a_k == b_k on either route ----
## An element that leaves the values of D and X unchanged commutes with the
## projector, so a_k == b_k for every b and p = 1 exactly: the no-power case
## of a layout whose d is constant within cells, or of a panel whose d and x
## vary by period only. .ipt_prepare() asserts that tie rather than
## recomputing it. The cell-level x below is in large units, which can send
## the projection to the QR route, whose rows for identical inputs may
## differ in the last bit; the tie must hold whichever route is taken, and
## row names on d and x must not hide it. (Without the check on the inputs
## the layout fit below gave p = 0.57 and a finite interval.)
set.seed(101)
n6 <- 6L; L6 <- 21L
row6 <- rep(rep(seq_len(n6), each = n6), each = L6)
col6 <- rep(rep(seq_len(n6), times = n6), each = L6)
cid6 <- (row6 - 1L) * n6 + col6
dc6 <- rnorm(n6 * n6); xc6 <- rnorm(n6 * n6)
y6 <- 0.3 * dc6[cid6] + 0.5 * xc6[cid6] + rnorm(length(cid6))
## one element by hand: reverse the replicates of every cell
idx6 <- split(seq_along(cid6), cid6)
g6 <- integer(length(cid6))
g6[unlist(idx6)] <- unlist(lapply(idx6, rev))
D6 <- matrix(dc6[cid6]); X6 <- cbind(1, xc6[cid6] * 1e5)
for (rn in list(NULL, paste0("r", seq_along(cid6)))) {
  rownames(D6) <- rn; rownames(X6) <- rn
  p6 <- mwperm:::.ipt_prepare(y6, D6, X6, list(seq_along(g6), g6))
  stopifnot(identical(p6$v, p6$u), identical(p6$W, p6$M))
}
## through the front ends: a layout, and a panel with a calendar-year column
f6 <- suppressWarnings(mwperm_layout(y6, dc6[cid6], xc6[cid6] * 1e5,
                                     row = row6, col = col6, seed = 1))
stopifnot(f6$pvalue == 1, identical(f6$conf_int, c(-Inf, Inf)))
set.seed(9001)
pd6 <- expand.grid(row = 1:21, col = 1:21, time = 1:8)
dt6 <- rnorm(8)
pd6$y <- 0.3 * dt6[pd6$time] + 0.05 * pd6$time + rnorm(nrow(pd6))
f6p <- suppressWarnings(mwperm_panel(pd6$y, dt6[pd6$time],
                                     x = 2000 + pd6$time, row = pd6$row,
                                     col = pd6$col, time = pd6$time,
                                     time_fe = FALSE, seed = 1))
stopifnot(f6p$pvalue == 1, identical(f6p$conf_int, c(-Inf, Inf)))

passed("test-projection.R")
