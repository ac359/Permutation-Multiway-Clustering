## Corrected reference implementation of the sign-flip test (IPT-Het).
##
## Sourced by tests/test-signflip.R, never run as a test (this directory is
## outside what R CMD check executes; see tests/README.md).
##
## This is a port of RPT_signflip() from the paper author's research script
## (utils.R, "additional functions in the revision"), kept as close to the
## original as possible -- the same array layout, the same explicit
## orthonormal complement through a complete QR, the same a/b statistics --
## so that agreement with the package certifies the PROCEDURE and not merely
## the package against itself. It differs from the original in exactly two
## authorized ways, both marked DEVIATION below and both confirmed with the
## author:
##
##   1. It enumerates one sign vector per coset {s, -s} (coordinate 1 pinned
##      to +1) instead of all 2^n_flip. Because the sign of observation (i, j)
##      is the PRODUCT s[g1(i)] * s[g2(j)], s and -s induce the identical
##      transformation, so the original computed every distinct element twice
##      -- every a_k and b_k appeared twice, the count and the group order both
##      doubled, and the p-value was unchanged. Deduplicating is therefore
##      exact, and halves the cost.
##
##   2. It excludes the identity element from the loop and from min_j a_j.
##      Equation (10) of Guo, Toulis and Wang (2026) takes the minimum over
##      1 <= j <= K, the NON-identity elements, and the Appendix D proof
##      defines the minorized function over that same range. The original
##      included the identity, where cbind(X, X_flip) has rank p rather than
##      2p and get_orthogonal_matrix() takes columns (2p+1):N of the complete
##      Q regardless -- a proper subspace of the true complement of col(X) --
##      so its a_0 could fall below every a_k and lower the minimum. Measured
##      on the author's setup (n = 12, n_flip = 5): min(a) lowered in ~7.5% of
##      draws, the p-value changed in ~2.5%, always in the conservative
##      direction, so the original was valid but was not Procedure 1.
##
## One further difference of convention, not of substance: the flip-group
## assignments g1 and g2 are ARGUMENTS here rather than drawn inside with
## sample(), so that the package and this reference can be handed the same
## assignment and compared on the same group. The package's own RNG hygiene
## (its assignment is seeded and leaves the caller's stream untouched) is
## tested separately.
##
## Arguments follow the original's array layout: Y is m x n x 1, X is
## m x n x 1 x p, D is m x n x 1. Flattening is column-major, as in the
## original's array2matrix(): observation (i, j) sits at position i + (j-1)*m.

## column-major flattening of a 3-d or 4-d array (the original's helper)
ref_array2matrix <- function(A) {
  if (length(dim(A)) == 4L)
    matrix(A, nrow = prod(dim(A)[1:3]), ncol = dim(A)[4])
  else matrix(A, nrow = prod(dim(A)[1:3]), ncol = 1L)
}

## orthonormal basis of the complement of col(X) via a complete QR (the
## original's helper; assumes X has full column rank, which every NON-identity
## stacked design [X | S_k X] has generically -- the identity's does not, which
## is deviation 2)
ref_orthogonal_matrix <- function(X) {
  Qfull <- qr.Q(qr(X), complete = TRUE)
  Qfull[, (ncol(X) + 1L):ncol(Qfull), drop = FALSE]
}

## The corrected reference. Returns the p-value for H0: beta = beta_null.
rpt_signflip_ref <- function(Y, X, D, group1, group2, n_flip,
                             beta_null = 0) {
  m <- dim(X)[1L]
  n <- dim(X)[2L]
  stopifnot(length(group1) == m, length(group2) == n,
            length(unique(c(group1, group2))) == n_flip)   # every group used

  ## DEVIATION 1: one representative per coset {s, -s}. The original built
  ##   R <- t(as.matrix(expand.grid(rep(list(c(-1, 1)), n_flip))))
  ## (n_flip x 2^n_flip) and looped over all 2^n_flip columns.
  R <- t(cbind(1, as.matrix(expand.grid(rep(list(c(-1, 1)), n_flip - 1L)))))
  n_elem <- ncol(R)                            # 2^(n_flip - 1)
  flip1 <- R[group1, , drop = FALSE]           # m x n_elem row signs
  flip2 <- R[group2, , drop = FALSE]           # n x n_elem column signs

  ## DEVIATION 2: the identity (every sign +1) is left out of the loop and of
  ## min(a). The original looped over it and included its a in the minimum.
  is_identity <- colSums(R == 1) == n_flip
  stopifnot(sum(is_identity) == 1L)
  ks <- which(!is_identity)

  X_mat <- ref_array2matrix(X)
  D_mat <- ref_array2matrix(D)
  ## Shift the outcome to the null being tested: under H0: beta = b the
  ## residual y - D b is what has the invariance. The original tested b = 0
  ## only; this is the only extension, and at b = 0 it is a no-op.
  Y0 <- Y - beta_null * D
  Y_mat <- ref_array2matrix(Y0)

  a <- numeric(length(ks))
  b <- numeric(length(ks))
  for (t in seq_along(ks)) {
    k <- ks[t]
    s1 <- flip1[, k]
    s2 <- flip2[, k]
    X_flip <- sweep(X,  MARGIN = 1L, STATS = s1, FUN = "*")
    X_flip <- sweep(X_flip, MARGIN = 2L, STATS = s2, FUN = "*")
    Y_flip <- sweep(Y0, MARGIN = 1L, STATS = s1, FUN = "*")
    Y_flip <- sweep(Y_flip, MARGIN = 2L, STATS = s2, FUN = "*")
    Xnew_mat <- cbind(X_mat, ref_array2matrix(X_flip))
    Vk <- ref_orthogonal_matrix(Xnew_mat)
    a[t] <- abs(t(D_mat) %*% Vk %*% t(Vk) %*% Y_mat)
    b[t] <- abs(t(D_mat) %*% Vk %*% t(Vk) %*% ref_array2matrix(Y_flip))
  }
  ## Equation (10): (1 + #{k : min_j a_j <= b_k}) / |G|, with |G| =
  ## 2^(n_flip - 1) the number of distinct elements INCLUDING the identity
  ## (the "1 +" is the identity's own contribution).
  (1 + sum(min(a) <= b)) / n_elem
}

## The heteroskedastic gravity DGP of the author's simulate_gravity_model(),
## in long format (one row per (i, j), i varying fastest -- the same order as
## the column-major flattening above). rho1 scales the error sd with the
## standardised gravity mean, rho2 with the standardised distance covariate;
## rho1 = rho2 = 0 is homoskedastic. The errors are centred lognormal
## (asymmetric), exactly as in the original; see the design note for why
## that matters for what the size check can and cannot show.
sim_gravity_het <- function(n, b, rho1, rho2, phi1 = 0.4, phi2 = 0.4) {
  m <- n
  log_gdp <- rnorm(n, mean = log(500), sd = 1.2)
  gdp <- pmin(pmax(exp(log_gdp), 5), 25000)
  log_gdp <- log(gdp)
  ## two-way random-effects covariate ("normal" cov_type of the original)
  s1sq <- phi1 / ((1 - phi1) * (1 - phi2) - phi1 * phi2)
  s2sq <- phi2 / ((1 - phi1) * (1 - phi2) - phi1 * phi2)
  vg <- rnorm(m); vh <- rnorm(n)
  g <- expand.grid(i = seq_len(m), j = seq_len(n))
  dist_raw <- sqrt(s1sq) * vg[g$i] + sqrt(s2sq) * vh[g$j] + rnorm(m * n)
  log_dist <- log(3000) + 0.5 * as.numeric(scale(dist_raw))
  X <- cbind(1, log_gdp[g$i], log_gdp[g$j])
  D <- log_dist
  mu <- 4 + 0.8 * X[, 2L] + 0.8 * X[, 3L] + b * D
  sigma <- 0.35 * exp(rho1 * as.numeric(scale(mu)) + rho2 * as.numeric(scale(D)))
  u <- rlnorm(m * n, meanlog = 0, sdlog = 1)
  eta <- (u - exp(0.5)) / sqrt((exp(1) - 1) * exp(1))   # centred, unit sd
  y <- mu + sigma * eta
  list(y = y, d = D, x = X[, 2:3], i = g$i, j = g$j, m = m, n = n)
}
