## Shared internals for all mwperm_* front ends. Not exported.

#' OLS reference estimate of the coefficient(s) of interest and a naive SE
#' (used only as a centre/scale for the confidence-interval search).
#' @keywords internal
#' @noRd
.ols_reference <- function(y, D, X) {
  D <- as.matrix(D); X <- as.matrix(X)
  d <- ncol(D); p <- ncol(X)
  W <- cbind(X, D)
  fit <- stats::lm.fit(W, y)
  cf <- fit$coefficients
  idx <- (p + 1L):(p + d)
  est <- cf[idx]
  ## naive homoskedastic SE
  res <- fit$residuals
  dfres <- length(y) - fit$rank
  se <- rep(NA_real_, d)
  if (dfres > 0L) {
    sigma2 <- sum(res^2) / dfres
    XtXi <- tryCatch(solve(crossprod(W)), error = function(e) NULL)
    if (!is.null(XtXi)) se <- sqrt(sigma2 * diag(XtXi)[idx])
  }
  list(estimate = as.numeric(est), se = as.numeric(se))
}

#' The IPT engine.
#'
#' @param perm_builder function(rep_seed) -> list of K+1 observation gather
#'   vectors (element 1 the identity), independent of any shift of `y`.
#' @keywords internal
#' @noRd
.ipt_engine <- function(y, D, X, perm_builder, K, n_reps, seed,
                        alpha, conf_int, beta_null, grid, conf_level,
                        type, d_names, n_clusters, call) {
  y <- as.numeric(y); D <- as.matrix(D); X <- as.matrix(X)
  N <- length(y); d <- ncol(D); p <- ncol(X)
  if (N <= 2L * p) {
    stop(sprintf(paste0("Need N > 2p for the projection to exist: N = %d, ",
                        "p = %d nuisance columns (incl. intercept). Drop ",
                        "covariates or add data."), N, p), call. = FALSE)
  }

  ## Reference OLS estimate / scale
  ref <- .ols_reference(y, D, X)

  ## The permuted-D projections (W) are needed for CI inversion / joint region
  ## inversion and, in the point p-value, whenever the null is non-zero (the
  ## b-statistic slope). They can be skipped only for a no-CI test of beta = 0
  ## (e.g. Monte-Carlo size).
  res_min <- 1 / (K + 1L)
  want_ci     <- isTRUE(conf_int) && d == 1L && res_min <= alpha
  want_region <- isTRUE(conf_int) && d >  1L && res_min <= alpha
  beta0 <- rep(beta_null, length.out = d)
  need_W <- want_ci || want_region || any(beta0 != 0)

  ## Build + prepare permutations once per rep (reused for the CI search)
  seeds <- if (is.null(seed)) rep(list(NULL), n_reps) else as.list(seed + seq_len(n_reps) - 1L)
  prep_list <- vector("list", n_reps)
  pv <- numeric(n_reps)
  for (r in seq_len(n_reps)) {
    op <- perm_builder(seeds[[r]])
    prep_list[[r]] <- .ipt_prepare(y, D, X, op, need_perm_D = need_W)
    pv[r] <- .ipt_eval(prep_list[[r]], beta0)$pvalue
  }
  pvalue <- stats::median(pv)
  Kp1 <- prep_list[[1L]]$Kp1

  ## Confidence set by test inversion: an interval for a single coefficient,
  ## a joint region (grid-based) for several.
  ci <- NULL; conf_region <- NULL; conf_box <- NULL
  note <- character(0)
  warned_res <- FALSE
  if (isTRUE(conf_int)) {
    if (res_min > alpha) {
      note <- c(note, sprintf(paste0("No %.0f%% confidence %s: the smallest attainable ",
                                     "p-value is 1/(K+1) = %.3g > alpha = %.3g, so every value ",
                                     "is retained. Increase K (more clusters/larger blocks)."),
                              100 * conf_level, if (d == 1L) "interval" else "region",
                              res_min, alpha))
      warned_res <- TRUE
    } else if (d == 1L) {
      ci <- .invert_ci(prep_list, alpha = 1 - conf_level,
                       centre = ref$estimate, scale = ref$se,
                       y = y, D = D, grid = grid)
    } else {
      reg <- .invert_region(prep_list, alpha = 1 - conf_level,
                            centre = ref$estimate, scale = ref$se,
                            d_names = d_names, grid = grid)
      conf_region <- reg$points; conf_box <- reg$box
      if (length(reg$note)) note <- c(note, reg$note)
    }
  }

  if (res_min > alpha && !warned_res) {
    note <- c(note, sprintf(paste0("Smallest attainable p-value is 1/(K+1) = %.3g > alpha = %.3g; ",
                                   "the test cannot reject at this level. Increase K (needs more ",
                                   "clusters per dimension)."), res_min, alpha))
  }

  structure(
    list(
      pvalue      = pvalue,
      pvalues_rep = pv,
      estimate    = stats::setNames(ref$estimate, d_names),
      se_naive    = stats::setNames(ref$se, d_names),
      conf_int    = ci,
      conf_region = conf_region,
      conf_box    = conf_box,
      conf_level  = conf_level,
      alpha       = alpha,
      beta_null   = beta_null,
      K           = Kp1 - 1L,
      n_perm      = Kp1,
      n_reps      = n_reps,
      type        = type,
      d_names     = d_names,
      n_obs       = N,
      n_clusters  = n_clusters,
      resolution  = res_min,
      note        = note,
      call        = call
    ),
    class = "mwperm"
  )
}

#' Test-inversion confidence interval for a single coefficient.
#'
#' Permutations are held fixed across candidate values of `b`, so the p-value
#' is a deterministic step function of `b` within each rep; the acceptance set
#' is taken to be an interval and its end points are located by outward
#' bracketing then bisection. The reported interval is the median over reps of
#' the per-rep interval end points. The p-value at each candidate `b` is read
#' off the cached prep objects (\code{\link{.ipt_prepare}}) in O(K), so the
#' whole search costs no extra QR decompositions.
#' @param prep_list list of per-rep prep objects (each with `has_perm_D = TRUE`).
#' @param y,D the (unshifted) outcome and single covariate, used only to derive
#'   a sensible step size when the naive SE is unavailable.
#' @keywords internal
#' @noRd
.invert_ci <- function(prep_list, alpha, centre, scale, y, D, grid = NULL,
                       max_expand = 60L, tol_factor = 1e-3) {
  if (!is.finite(centre)) centre <- 0
  step <- scale
  if (!is.finite(step) || step <= 0) {
    step <- stats::sd(y) / max(stats::sd(as.numeric(D)), .Machine$double.eps)
    if (!is.finite(step) || step <= 0) step <- 1
  }
  tol <- step * tol_factor

  pval_at <- function(b, prep) .ipt_eval(prep, b)$pvalue

  one_side <- function(prep, direction) {
    ## centre (the point estimate) rejected => degenerate; fall back to centre
    if (pval_at(centre, prep) <= alpha) return(centre)
    lo <- centre                       # last value known to be accepted
    hi <- NA_real_                     # first value known to be rejected
    h <- step
    for (i in seq_len(max_expand)) {
      cand <- centre + direction * h
      if (pval_at(cand, prep) <= alpha) { hi <- cand; break }
      lo <- cand                       # cand accepted: advance the bracket
      h <- h * 1.6
    }
    if (is.na(hi)) return(direction * Inf)  # never rejected: unbounded this side
    ## bisect, maintaining pval(lo) > alpha >= pval(hi)
    while (abs(hi - lo) > tol) {
      mid <- (lo + hi) / 2
      if (pval_at(mid, prep) > alpha) lo <- mid else hi <- mid
    }
    lo
  }

  if (!is.null(grid)) {
    ## explicit grid mode: union over reps of accepted grid points
    acc <- rep(FALSE, length(grid))
    for (prep in prep_list) {
      pg <- vapply(grid, pval_at, numeric(1), prep = prep)
      acc <- acc | (pg > alpha)
    }
    if (!any(acc)) return(c(NA_real_, NA_real_))
    return(range(grid[acc]))
  }

  lowers <- numeric(length(prep_list)); uppers <- numeric(length(prep_list))
  for (r in seq_along(prep_list)) {
    lowers[r] <- one_side(prep_list[[r]], -1)
    uppers[r] <- one_side(prep_list[[r]], +1)
  }
  c(stats::median(lowers), stats::median(uppers))
}

#' Joint confidence region for several coefficients by test inversion.
#'
#' Inverts the exact joint test \eqn{H_0: \beta = b} over a grid of candidate
#' vectors `b`: each retained point is one at which the (median-over-reps) test
#' does not reject at level `alpha`, so the retained set is a finite-sample
#' valid \eqn{1 - \alpha} confidence region (discretised by the grid). Cheap
#' because every evaluation reuses the cached `prep` objects via
#' \code{\link{.ipt_eval}} (no QR refactorisation).
#'
#' @param prep_list list of per-rep prep objects (with `has_perm_D = TRUE`).
#' @param centre,scale length-d OLS estimate / naive SE, used to place a default
#'   grid around the estimate.
#' @param d_names column labels for the returned box.
#' @param grid optional explicit grid: a list of d numeric vectors (one set of
#'   candidate values per coordinate), or a single vector used for every
#'   coordinate. When `NULL`, a default grid of `n_grid` points spanning
#'   `centre +/- spread * scale` per coordinate is used.
#' @return list(points, box, note): `points` is the matrix of accepted `beta`
#'   vectors, `box` a 2 x d matrix of marginal (lower, upper) extents.
#' @keywords internal
#' @noRd
.invert_region <- function(prep_list, alpha, centre, scale, d_names, grid = NULL,
                           n_grid = 21L, spread = 6, max_points = 2e4L,
                           max_expand = 6L) {
  d <- length(centre)
  ck <- ifelse(is.finite(centre), centre, 0)
  sk <- ifelse(is.finite(scale) & scale > 0, scale, 1)

  ## evaluate acceptance over a Cartesian grid (median p-value over reps > alpha)
  eval_grid <- function(axes) {
    if (prod(vapply(axes, length, numeric(1))) > max_points) return(NULL)  # too big
    G <- as.matrix(expand.grid(axes))
    acc <- vapply(seq_len(nrow(G)), function(i) {
      b <- G[i, ]
      stats::median(vapply(prep_list, function(pp) .ipt_eval(pp, b)$pvalue,
                           numeric(1))) > alpha
    }, logical(1))
    list(G = G, acc = acc, axes = axes)
  }
  touches <- function(res) {
    if (!any(res$acc)) return(FALSE)
    pts <- res$G[res$acc, , drop = FALSE]
    any(vapply(seq_len(d), function(k)
      min(pts[, k]) <= min(res$axes[[k]]) || max(pts[, k]) >= max(res$axes[[k]]),
      logical(1)))
  }

  if (is.null(grid)) {
    ## default grid, expanded outward until the region no longer clips the edge
    sp <- spread
    res <- eval_grid(lapply(seq_len(d), function(k)
      seq(ck[k] - sp * sk[k], ck[k] + sp * sk[k], length.out = n_grid)))
    if (is.null(res))
      return(list(points = NULL, box = NULL, note = sprintf(paste0(
        "Joint confidence region skipped: a %d^%d default grid is too large. ",
        "Pass an explicit `grid` (a list of per-coefficient value vectors)."),
        n_grid, d)))
    for (e in seq_len(max_expand)) {
      if (!any(res$acc) || !touches(res)) break
      sp <- sp * 1.8
      r2 <- eval_grid(lapply(seq_len(d), function(k)
        seq(ck[k] - sp * sk[k], ck[k] + sp * sk[k], length.out = n_grid)))
      if (is.null(r2)) break
      res <- r2
    }
  } else {
    axes <- if (is.list(grid)) {
      if (length(grid) != d)
        stop(sprintf("`grid` list must have one vector per coefficient (d = %d).", d),
             call. = FALSE)
      grid
    } else rep(list(grid), d)
    res <- eval_grid(axes)
    if (is.null(res))
      return(list(points = NULL, box = NULL, note = sprintf(paste0(
        "Joint confidence region skipped: `grid` has > %d points."), max_points)))
  }

  if (!any(res$acc))
    return(list(points = res$G[0, , drop = FALSE], box = NULL,
                note = "Joint confidence region is empty on the searched grid; pass a finer/shifted `grid`."))

  pts <- res$G[res$acc, , drop = FALSE]
  colnames(pts) <- d_names
  box <- rbind(apply(pts, 2L, min), apply(pts, 2L, max))
  dimnames(box) <- list(c("lower", "upper"), d_names)
  note <- if (touches(res))
    "Joint confidence region still reaches the grid boundary (it may be unbounded); pass an explicit `grid` to widen it."
  else character(0)
  list(points = pts, box = box, note = note)
}

#' Validate and coerce a cluster id vector to dense 1-based integers.
#' @keywords internal
#' @noRd
.dense_id <- function(x) {
  f <- as.integer(factor(x))
  f
}

#' Error if any supplied vector does not have length N.
#'
#' @param N expected length (the number of observations).
#' @param vars a named list of vectors to check; names appear in the message.
#' @keywords internal
#' @noRd
.check_lengths <- function(N, vars) {
  bad <- names(vars)[vapply(vars, length, integer(1L)) != N]
  if (length(bad))
    stop(sprintf("%s must have the same length as `y` (%d).",
                 paste0("`", bad, "`", collapse = ", "), N), call. = FALSE)
  invisible(NULL)
}

#' Assemble the nuisance design X (with intercept) from a covariate spec.
#' @keywords internal
#' @noRd
.make_X <- function(x, N, intercept = TRUE) {
  if (is.null(x)) {
    X <- matrix(numeric(0), nrow = N, ncol = 0)
  } else {
    X <- as.matrix(x)
    if (nrow(X) != N) stop("`x` must have the same number of rows as `y`.", call. = FALSE)
  }
  if (intercept) X <- cbind(`(Intercept)` = 1, X)
  X
}
