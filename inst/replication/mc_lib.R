## inst/replication/mc_lib.R -- shared Monte Carlo harness for the mwperm
## replication scripts. Sourced by every 0X_*.R script in this directory.
##
## Provides:
##   * rf_array()        -- the random-feature error/covariate model
##                          (Guo, Toulis & Wang 2026, Eq. 13)
##   * dgp_dyadic71()    -- the paper's Section 7.1 dyadic DGP
##   * dgp_panel(), dgp_threeway() -- panel (arbitrary time trend, condition
##                          InvB) and fully exchangeable three-way (InvA) DGPs
##   * naive_ols_p()     -- classical OLS t-test p-value (the invalid baseline
##                          the permutation test is meant to replace)
##   * mc_cell()         -- batched, checkpointed, resumable, parallel runner
##   * size_table(), superuniformity(), fixed_u_check(), cp_* -- reporting
##                          helpers with Monte-Carlo SEs and Clopper-Pearson
##                          bounds
##
## Reproducibility conventions (identical to the package's own test scheme):
##   - Seeds: sim s in a cell draws its data at DGP seed base_seed + s and
##     fits at seed 100000 + s. Results are deterministic given (cell_id,
##     base_seed); parallel scheduling cannot change them because each sim
##     reseeds itself. Fit seeds must stay below 2,147,483 (the package's
##     .sub_seed = seed*1000 + j overflows integers above that).
##   - Cache: <cache_root>/<cell_id>/batch_%04d.rds + params.rds, where
##     cache_root defaults to "cache" under the working directory and can be
##     redirected with the MWPERM_REPL_CACHE environment variable. A params
##     mismatch on an existing cell is a hard error (never a silent reuse):
##     bump the cell_id's version tag or delete the cache directory.

MC_CACHE_ROOT <- Sys.getenv("MWPERM_REPL_CACHE", "cache")
MC_CORES <- max(1L, min(8L, parallel::detectCores() - 1L))

## ---- Eq. (13): random-feature array ----------------------------------------
## variable_ij = sigma1 * v_{1,i} + sigma2 * v_{2,j} + v_{3,ij}, with
## sigma1^2 = phi1/(1 - phi1 - phi2), sigma2^2 = phi2/(1 - phi1 - phi2),
## inducing correlation phi1 along dimension 1 and phi2 along dimension 2.
rf_array <- function(n1, n2, phi1, phi2, rdist = stats::rnorm) {
  stopifnot(phi1 >= 0, phi2 >= 0, phi1 + phi2 < 1)
  s1 <- sqrt(phi1 / (1 - phi1 - phi2))
  s2 <- sqrt(phi2 / (1 - phi1 - phi2))
  outer(s1 * rdist(n1), s2 * rdist(n2), "+") + matrix(rdist(n1 * n2), n1, n2)
}

## ---- Section 7.1 dyadic DGP -------------------------------------------------
## p = 3, x_ij = (1, z_i, z_j), z ~ Unif[0, 2], gamma = (0.5, 1, 1); the
## covariate w and the error come from Eq. (13); d = w (normal) or
## d = exp(0.5 w) (lognormal). Column-major stacking (i varies fastest).
dgp_dyadic71 <- function(n, cov_type = c("normal", "lognormal"),
                         phi2_err, beta = 0) {
  cov_type <- match.arg(cov_type)
  z  <- stats::runif(n, 0, 2)
  ii <- rep(seq_len(n), times = n)
  jj <- rep(seq_len(n), each  = n)
  w  <- rf_array(n, n, 0.4, 0.4)
  d  <- if (cov_type == "lognormal") exp(0.5 * as.vector(w)) else as.vector(w)
  eps <- as.vector(rf_array(n, n, 0.05, phi2_err))
  x  <- cbind(z_i = z[ii], z_j = z[jj])          # front end adds the intercept
  y  <- 0.5 + x %*% c(1, 1) + d * beta + eps
  list(y = as.vector(y), d = d, x = x, row = ii, col = jj)
}

## ---- panel DGP (condition InvB): arbitrary common time trend ----------------
## y_ijt = 1.5 + x'gamma + d*beta + eta_i + xi_j + zeta_t + u_ijt, with zeta_t
## an arbitrary (irregular) trend shared across all pairs -- errors are
## exchangeable across (i, j) within each period but NOT over time. A random
## treatment-timing covariate d supplies the coefficient of interest.
dgp_panel <- function(n, Tt, beta = 0, trend = TRUE) {
  z  <- stats::runif(n, 0, 2)
  ii <- rep(rep(seq_len(n), times = n), times = Tt)
  jj <- rep(rep(seq_len(n), each  = n), times = Tt)
  tt <- rep(seq_len(Tt), each = n * n)
  eta <- stats::rnorm(n); xi <- stats::rnorm(n)
  zeta <- if (trend) cumsum(stats::rnorm(Tt, sd = 1.5)) else stats::rnorm(Tt)
  ## treatment turns on at a pair-specific random date (0/1, dyad-by-time)
  on <- matrix(stats::runif(n * n) < 0.5, n, n)
  start <- matrix(sample.int(Tt, n * n, replace = TRUE), n, n)
  d <- as.numeric(vapply(seq_along(tt), function(k)
    on[ii[k], jj[k]] && tt[k] >= start[ii[k], jj[k]], logical(1)))
  u <- stats::rnorm(n * n * Tt)
  x <- cbind(z_i = z[ii], z_j = z[jj])
  y <- 1.5 + x %*% c(0.6, 0.6) + d * beta +
    eta[ii] + xi[jj] + zeta[tt] + u
  list(y = as.vector(y), d = d, x = x, row = ii, col = jj, time = tt)
}

## ---- three-way DGP (condition InvA): fully exchangeable, no trend -----------
## y_ijt = 1 + x'gamma + d*beta + eta_i + xi_j + rho_t + u_ijt, every dimension
## exchangeable (rho_t i.i.d., no ordering) -- the setting mwperm_threeway is
## valid for. Feeding dgp_panel(trend = TRUE) to the three-way test instead is
## the negative control in 04_negative_control.R.
dgp_threeway <- function(n, Tt, beta = 0) {
  z  <- stats::runif(n, 0, 2)
  ii <- rep(rep(seq_len(n), times = n), times = Tt)
  jj <- rep(rep(seq_len(n), each  = n), times = Tt)
  tt <- rep(seq_len(Tt), each = n * n)
  eta <- stats::rnorm(n); xi <- stats::rnorm(n); rho <- stats::rnorm(Tt)
  d <- stats::rnorm(n * n * Tt)
  u <- stats::rnorm(n * n * Tt)
  x <- cbind(z_i = z[ii], z_j = z[jj])
  y <- 1 + x %*% c(0.6, 0.6) + d * beta + eta[ii] + xi[jj] + rho[tt] + u
  list(y = as.vector(y), d = d, x = x, id1 = ii, id2 = jj, id3 = tt)
}

## ---- classical OLS t-test (the invalid baseline) ----------------------------
## Two-sided p-value for the last column of the design under i.i.d.-error OLS
## inference. On multi-way clustered data this over-rejects badly -- it is the
## failure mode the permutation test is designed to avoid, shown side by side.
naive_ols_p <- function(y, X, d) {
  W <- cbind(X, d)                     # X already includes an intercept column
  fit <- stats::lm.fit(W, y)
  k <- ncol(W)
  dfres <- length(y) - fit$rank
  sigma2 <- sum(fit$residuals^2) / dfres
  XtXi <- tryCatch(chol2inv(chol(crossprod(W))), error = function(e) NULL)
  if (is.null(XtXi)) return(NA_real_)
  se <- sqrt(sigma2 * diag(XtXi)[k])
  unname(2 * stats::pt(-abs(fit$coefficients[k] / se), dfres))
}

## ---- batched, checkpointed Monte Carlo runner ------------------------------
## sim_fun(s, dgp_seed, fit_seed) must return a named numeric vector of fixed
## length. It runs inside a forked worker; it reseeds itself with dgp_seed
## before drawing data and passes fit_seed to the mwperm call.
mc_cell <- function(cell_id, n_sims, sim_fun, params = list(),
                    batch = 500L, base_seed = 20260713L,
                    cores = MC_CORES, cache_root = MC_CACHE_ROOT,
                    quiet = FALSE) {
  dir <- file.path(cache_root, cell_id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  pfile <- file.path(dir, "params.rds")
  params_full <- c(params, list(.base_seed = base_seed, .batch = batch))
  if (file.exists(pfile)) {
    if (!identical(readRDS(pfile), params_full))
      stop("cell '", cell_id, "': cached params.rds does not match current ",
           "parameters. Bump the cell_id version tag or delete the cache dir.")
  } else saveRDS(params_full, pfile)

  n_batches <- ceiling(n_sims / batch)
  out <- vector("list", n_batches)
  for (b in seq_len(n_batches)) {
    f <- file.path(dir, sprintf("batch_%04d.rds", b))
    idx <- seq.int((b - 1L) * batch + 1L, min(b * batch, n_sims))
    if (file.exists(f)) {
      m <- readRDS(f)
      if (nrow(m) >= length(idx)) {
        out[[b]] <- m[seq_along(idx), , drop = FALSE]
        next
      }
    }
    t0 <- Sys.time()
    res <- parallel::mclapply(idx, function(s) {
      tryCatch({
        set.seed(base_seed + s)
        sim_fun(s, dgp_seed = base_seed + s, fit_seed = 100000L + s)
      }, error = function(e)
        structure(list(sim = s, msg = conditionMessage(e)),
                  class = "mc_error"))
    }, mc.cores = cores, mc.preschedule = TRUE)
    bad <- vapply(res, inherits, logical(1), what = "mc_error")
    if (any(bad)) {
      err <- res[bad][[1]]
      stop("cell '", cell_id, "', sim ", err$sim, " failed: ", err$msg)
    }
    m <- do.call(rbind, res)
    attr(m, "elapsed_s") <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    saveRDS(m, f)
    out[[b]] <- m
    if (!quiet)
      cat(sprintf("  [%s] batch %d/%d (%d sims) in %.1fs\n",
                  cell_id, b, n_batches, length(idx), attr(m, "elapsed_s")))
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res[seq_len(n_sims), , drop = FALSE]
}

## ---- Clopper-Pearson helpers -----------------------------------------------
cp_interval <- function(x, n, conf = 0.95) {
  a <- 1 - conf
  c(lo = if (x == 0) 0 else stats::qbeta(a / 2, x, n - x + 1),
    hi = if (x == n) 1 else stats::qbeta(1 - a / 2, x + 1, n - x))
}
cp_upper1 <- function(x, n, conf = 0.95)
  if (x == n) 1 else stats::qbeta(conf, x + 1, n - x)
cp_lower1 <- function(x, n, conf = 0.95)
  if (x == 0) 0 else stats::qbeta(1 - conf, x, n - x + 1)

## ---- size table -------------------------------------------------------------
## Rejection rule is p <= alpha (the package's convention). Reports the
## estimate, its Monte-Carlo SE, the two-sided 95% CP interval, and the
## one-sided upper bound used to certify size control.
size_table <- function(p, alphas = c(0.01, 0.05, 0.10)) {
  p <- p[!is.na(p)]
  n <- length(p)
  do.call(rbind, lapply(alphas, function(a) {
    x <- sum(p <= a)
    ci <- cp_interval(x, n)
    data.frame(alpha = a, n_sims = n, n_reject = x, size = x / n,
               mc_se = sqrt((x / n) * (1 - x / n) / n),
               cp_lo = ci["lo"], cp_hi = ci["hi"],
               cp_upper1 = cp_upper1(x, n), row.names = NULL)
  }))
}

## ---- super-uniformity over the attainable p-grid ---------------------------
## Exactness means P(p <= u) <= u for every u; p lives on {1,...,K+1}/(K+1).
## For each grid point: the ECDF and a Bonferroni-corrected one-sided CP lower
## bound; a point whose corrected lower bound exceeds u is a validity
## violation.
superuniformity <- function(p, K, conf_family = 0.95) {
  p <- p[!is.na(p)]
  stopifnot(length(unique(K)) == 1)
  Kp1 <- unique(K) + 1L
  u <- seq_len(Kp1) / Kp1
  n <- length(p)
  conf_pt <- 1 - (1 - conf_family) / Kp1
  out <- do.call(rbind, lapply(u, function(uu) {
    x <- sum(p <= uu + 1e-12)
    data.frame(u = uu, ecdf = x / n, excess = x / n - uu,
               cp_lower1 = cp_lower1(x, n, conf_pt),
               viol = cp_lower1(x, n, conf_pt) > uu + 1e-12, row.names = NULL)
  }))
  attr(out, "any_violation") <- any(out$viol)
  out
}

## ---- fixed-u validity check (for cells where K varies across sims) ----------
fixed_u_check <- function(p, u_set = c(0.05, 0.10, 0.20, 0.50),
                          conf_family = 0.95) {
  p <- p[!is.na(p)]
  n <- length(p)
  conf_pt <- 1 - (1 - conf_family) / length(u_set)
  out <- do.call(rbind, lapply(u_set, function(uu) {
    x <- sum(p <= uu + 1e-12)
    data.frame(u = uu, ecdf = x / n, excess = x / n - uu,
               cp_lower1 = cp_lower1(x, n, conf_pt),
               viol = cp_lower1(x, n, conf_pt) > uu + 1e-12, row.names = NULL)
  }))
  attr(out, "any_violation") <- any(out$viol)
  out
}

## ---- misc -------------------------------------------------------------------
fmt_pct <- function(x, d = 1) sprintf(paste0("%.", d, "f"), 100 * x)

sink_both <- function(file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  sink(file, split = TRUE)
  invisible(file)
}
