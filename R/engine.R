## Shared internals for all mwperm_* front ends. Not exported.

#' OLS reference estimate of the coefficient(s) of interest and a naive SE
#' (used only as a centre/scale for the confidence-interval search).
#' @keywords internal
#' @noRd
.ols_reference <- function(y, D, X) {
  D <- as.matrix(D); X <- as.matrix(X)
  d <- ncol(D)                         # number of coefficients of interest
  p <- ncol(X)                         # number of nuisance columns (incl. intercept)
  W <- cbind(X, D)                     # full design: [nuisance | covariate(s) of interest]
  fit <- stats::lm.fit(W, y)
  cf <- fit$coefficients
  idx <- (p + 1L):(p + d)              # positions of the D coefficients within `cf`
  est <- cf[idx]                       # OLS point estimate(s) of beta
  ## naive homoskedastic SE -- used only as a centre/scale for the CI search,
  ## never reported as an inferential quantity.
  res <- fit$residuals
  dfres <- length(y) - fit$rank        # residual degrees of freedom
  se <- rep(NA_real_, d)
  if (dfres > 0L) {
    sigma2 <- sum(res^2) / dfres       # homoskedastic error-variance estimate
    XtXi <- tryCatch(solve(crossprod(W)), error = function(e) NULL)  # (W'W)^{-1}, or NULL if singular
    if (!is.null(XtXi)) se <- sqrt(sigma2 * diag(XtXi)[idx])
  }
  list(estimate = as.numeric(est), se = as.numeric(se))
}

#' The IPT engine: shared driver behind every mwperm_* front end.
#'
#' Given the data and a design-specific permutation builder, this runs the test
#' for `n_reps` independent random permutation groups, aggregates the p-value by
#' the median, and (optionally) inverts the test to a confidence set. All five
#' front ends differ only in how they build their permutations and validate
#' their inputs; everything downstream of that is handled here.
#'
#' @param y,D,X numeric outcome, covariate(s) of interest, and nuisance design
#'   (intercept already included), all with N rows.
#' @param perm_builder function(rep_seed) -> list of K+1 observation gather
#'   vectors (element 1 the identity), independent of any shift of `y`.
#' @param K,n_reps,seed group order (non-identity count), number of repetitions,
#'   and base RNG seed.
#' @param alpha,conf_int,beta_null,grid test level, whether to invert a
#'   confidence set, the null value(s), and an optional inversion grid.
#' @param type,d_names,n_clusters,call metadata stored on the result for
#'   printing.
#' @return an object of class `"mwperm"`: a list whose fields are documented
#'   inline at the `structure()` call below and consumed by the S3 methods in
#'   methods.R.
#' @keywords internal
#' @noRd
.ipt_engine <- function(y, D, X, perm_builder, K, n_reps, seed,
                        alpha, conf_int, beta_null, grid,
                        type, d_names, n_clusters, call) {
  y <- as.numeric(y); D <- as.matrix(D); X <- as.matrix(X)
  N <- length(y)                     # number of observations
  d <- ncol(D)                       # number of coefficients of interest
  p <- ncol(X)                       # number of nuisance columns (incl. intercept)
  conf_level <- 1 - alpha            # always the complement of the test level
  ## Need more observations than the stacked projection [X | X_k] consumes; with
  ## p nuisance columns that projection has up to 2p columns.
  if (N <= 2L * p) {
    stop(sprintf(paste0("Need N > 2p for the projection to exist: N = %d, ",
                        "p = %d nuisance columns (incl. intercept). Drop ",
                        "covariates or add data."), N, p), call. = FALSE)
  }

  ## Reference OLS estimate / scale (centre + step size for the CI search)
  ref <- .ols_reference(y, D, X)

  ## res_min: the p-value resolution, i.e. the smallest attainable p-value
  ## 1/(K+1). If it already exceeds alpha no confidence set can exclude anything.
  res_min <- 1 / (K + 1L)
  ## The permuted-D projections (W in the prep object) are needed for CI / joint
  ## region inversion and, in the point p-value, whenever the null is non-zero
  ## (they carry the b-statistic slope). They can be skipped only for a no-CI
  ## test of beta = 0 (e.g. a Monte-Carlo size simulation) to save work.
  want_ci     <- isTRUE(conf_int) && d == 1L && res_min <= alpha   # scalar interval feasible
  want_region <- isTRUE(conf_int) && d >  1L && res_min <= alpha   # joint region feasible
  beta0 <- rep(beta_null, length.out = d)                          # null value(s), recycled to length d
  need_W <- want_ci || want_region || any(beta0 != 0)

  ## Build + prepare the permutations once per rep, then reuse the cached prep
  ## objects for the CI search (so the QR work is never repeated).
  ## seeds: one per rep (NULL = use the ambient RNG, no reproducibility).
  seeds <- if (is.null(seed)) rep(list(NULL), n_reps) else as.list(seed + seq_len(n_reps) - 1L)
  prep_list <- vector("list", n_reps)   # cached cross products, one per rep
  pv <- numeric(n_reps)                 # per-rep p-value at the null
  for (r in seq_len(n_reps)) {
    op <- perm_builder(seeds[[r]])      # K+1 observation gather-vectors for this rep
    prep_list[[r]] <- .ipt_prepare(y, D, X, op, need_perm_D = need_W)
    pv[r] <- .ipt_eval(prep_list[[r]], beta0)$pvalue
  }
  pvalue <- stats::median(pv)           # reported p-value: median across reps
  Kp1 <- prep_list[[1L]]$Kp1            # realised group order (K + 1)

  ## Confidence set by test inversion: an interval for a single coefficient,
  ## a joint region (grid-based) for several.
  ci <- NULL; conf_region <- NULL; conf_box <- NULL
  note <- character(0)                  # human-readable caveats appended below
  warned_res <- FALSE                   # guard so the coarse-resolution note is added once
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

  ## Assemble the returned "mwperm" object. Field meanings (read by the S3
  ## methods in methods.R):
  structure(
    list(
      pvalue      = pvalue,        # reported p-value (median over reps)
      pvalues_rep = pv,            # the per-rep p-values (plotted by plot.mwperm)
      estimate    = stats::setNames(ref$estimate, d_names),  # OLS point estimate(s)
      se_naive    = stats::setNames(ref$se, d_names),        # naive homoskedastic SE(s)
      conf_int    = ci,            # length-2 inverted interval (d == 1), else NULL
      conf_region = conf_region,   # matrix of retained beta vectors (d > 1), else NULL
      conf_box    = conf_box,      # 2 x d marginal extent of the region, else NULL
      conf_level  = conf_level,    # 1 - alpha
      alpha       = alpha,         # test level
      beta_null   = beta_null,     # null value tested
      K           = Kp1 - 1L,      # number of non-identity permutations
      n_perm      = Kp1,           # group order (K + 1)
      n_reps      = n_reps,        # number of independent repetitions aggregated
      type        = type,          # design label, e.g. "dyadic"
      d_names     = d_names,       # coefficient name(s)
      n_obs       = N,             # observations actually used
      n_clusters  = n_clusters,    # named per-dimension cluster counts
      resolution  = res_min,       # p-value resolution 1/(K+1)
      note        = note,          # character vector of caveats
      call        = call           # the originating front-end call
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
  ## centre: where the bracketing starts (the OLS point estimate).
  if (!is.finite(centre)) centre <- 0
  ## step: initial bracketing increment. Prefer the naive SE; if it is
  ## unavailable fall back to a crude scale sd(y)/sd(D), then to 1.
  step <- scale
  if (!is.finite(step) || step <= 0) {
    step <- stats::sd(y) / max(stats::sd(as.numeric(D)), .Machine$double.eps)
    if (!is.finite(step) || step <= 0) step <- 1
  }
  tol <- step * tol_factor             # bisection stopping width

  pval_at <- function(b, prep) .ipt_eval(prep, b)$pvalue   # p-value at null b, one rep

  ## Locate one endpoint of the acceptance interval for a single rep.
  ## `direction` is -1 (lower limit) or +1 (upper limit).
  one_side <- function(prep, direction) {
    ## centre (the point estimate) rejected => degenerate; fall back to centre
    if (pval_at(centre, prep) <= alpha) return(centre)
    lo <- centre                       # last value known to be accepted
    hi <- NA_real_                     # first value known to be rejected
    h <- step                          # current step out from the centre
    for (i in seq_len(max_expand)) {   # phase 1: expand outward to bracket the edge
      cand <- centre + direction * h
      if (pval_at(cand, prep) <= alpha) { hi <- cand; break }
      lo <- cand                       # cand accepted: advance the bracket
      h <- h * 1.6                     # geometric growth so few steps are needed
    }
    if (is.na(hi)) return(direction * Inf)  # never rejected: unbounded this side
    ## phase 2: bisect, maintaining pval(lo) > alpha >= pval(hi)
    while (abs(hi - lo) > tol) {
      mid <- (lo + hi) / 2
      if (pval_at(mid, prep) > alpha) lo <- mid else hi <- mid
    }
    lo                                 # the accepted side of the converged bracket
  }

  if (!is.null(grid)) {
    ## explicit grid mode: a point is retained if accepted in ANY rep (union),
    ## and the reported interval is the range of retained points.
    acc <- rep(FALSE, length(grid))    # accepted-anywhere flag per grid point
    for (prep in prep_list) {
      pg <- vapply(grid, pval_at, numeric(1), prep = prep)
      acc <- acc | (pg > alpha)
    }
    if (!any(acc)) return(c(NA_real_, NA_real_))
    return(range(grid[acc]))
  }

  ## default mode: bracket each side in every rep, report the median endpoints.
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
  d <- length(centre)                              # number of coefficients
  ck <- ifelse(is.finite(centre), centre, 0)       # per-coord grid centre (0 if estimate missing)
  sk <- ifelse(is.finite(scale) & scale > 0, scale, 1)  # per-coord grid scale (1 if SE missing)

  ## Evaluate acceptance over a Cartesian grid: a point b is retained when its
  ## median p-value across reps exceeds alpha. `axes` is a list of d coordinate
  ## value-vectors; returns NULL (caller treats as "too large") if the product
  ## grid would exceed max_points.
  eval_grid <- function(axes) {
    if (prod(vapply(axes, length, numeric(1))) > max_points) return(NULL)  # too big
    G <- as.matrix(expand.grid(axes))             # all candidate beta vectors (rows)
    acc <- vapply(seq_len(nrow(G)), function(i) {  # retained-flag per candidate
      b <- G[i, ]
      stats::median(vapply(prep_list, function(pp) .ipt_eval(pp, b)$pvalue,
                           numeric(1))) > alpha
    }, logical(1))
    list(G = G, acc = acc, axes = axes)
  }
  ## TRUE if any retained point sits on the outer edge of the grid, i.e. the
  ## region was clipped and the grid should be widened.
  touches <- function(res) {
    if (!any(res$acc)) return(FALSE)
    pts <- res$G[res$acc, , drop = FALSE]         # retained points only
    any(vapply(seq_len(d), function(k)
      min(pts[, k]) <= min(res$axes[[k]]) || max(pts[, k]) >= max(res$axes[[k]]),
      logical(1)))
  }

  if (is.null(grid)) {
    ## default grid: centre +/- sp*scale per coordinate, expanded outward (sp
    ## grows) until the retained region no longer clips the grid edge.
    sp <- spread                       # current half-width in scale units
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
    ## user-supplied grid: a list of one value-vector per coefficient, or a
    ## single vector reused for every coefficient.
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

  pts <- res$G[res$acc, , drop = FALSE]            # retained beta vectors
  colnames(pts) <- d_names
  ## box: per-coordinate (lower, upper) extent of the retained set
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

## ---- shared front-end helpers --------------------------------------------
## Used by every mwperm_* front end (dyadic, panel, threeway, layout, missing)
## to pick the permutation-group order, derive coefficient labels, derive
## per-rep seeds, and validate complete-array designs.

#' Default / validate the permutation-group order K.
#'
#' The group has order \code{K + 1}, so the smallest attainable p-value is
#' \code{1 / (K + 1)} and a non-trivial group needs \code{K + 1 <= } the
#' smallest permuted dimension. When \code{K} is \code{NULL} the largest group
#' the design supports is used, capped at \code{cap}.
#'
#' @param K user-supplied \code{K} or \code{NULL}.
#' @param dim_sizes the sizes of the permuted dimensions (scalar or vector); the
#'   smallest one bounds the group order.
#' @keywords internal
#' @noRd
.default_K <- function(K, dim_sizes, cap = 199L) {
  smallest <- min(dim_sizes)
  gmax <- smallest - 1L
  if (is.null(K)) {
    K <- min(gmax, cap)
  } else {
    K <- as.integer(K)
    if (K + 1L > smallest)
      stop(sprintf("K + 1 = %d exceeds the smallest permuted dimension (%d).",
                   K + 1L, smallest), call. = FALSE)
  }
  if (K < 1L) stop("Not enough clusters to permute (need >= 2 in each dimension).",
                   call. = FALSE)
  K
}

#' Derive a per-rep, per-dimension seed from a rep-level seed (or NULL).
#' @keywords internal
#' @noRd
.sub_seed <- function(rep_seed, j) if (is.null(rep_seed)) NULL else rep_seed * 1000L + j

#' Column labels for the coefficient(s) of interest.
#'
#' Uses the column names of \code{D} when present, else the deparsed user
#' expression \code{fallback} (suffixed by column index when \code{D} has
#' several columns).
#' @keywords internal
#' @noRd
.coef_names <- function(D, fallback) {
  nm <- colnames(D)
  if (!is.null(nm)) return(nm)
  if (ncol(D) == 1L) fallback else paste0(fallback, seq_len(ncol(D)))
}

#' Error unless a design is a complete balanced array (one observation per cell,
#' every cell present).
#'
#' Consolidates the identical check used by \code{\link{mwperm_panel}} and
#' \code{\link{mwperm_threeway}}.
#'
#' @param coords integer matrix of cluster coordinates, one row per observation.
#' @param sizes named integer vector of per-dimension sizes; the names label the
#'   dimensions in the error messages.
#' @param N the number of observations.
#' @param what a short noun phrase naming the design (for the message).
#' @keywords internal
#' @noRd
.require_complete_array <- function(coords, sizes, N, what) {
  dims <- paste(names(sizes), collapse = ", ")
  if (anyDuplicated(coords))
    stop(sprintf("Each (%s) cell must appear at most once.", dims), call. = FALSE)
  expected <- prod(sizes)
  if (N != expected)
    stop(sprintf(paste0("%s must be a complete balanced array: expected %d ",
                        "cells (%s) but found %d. Fill or drop cells so the ",
                        "array is complete."),
                 what, expected,
                 paste(sprintf("%s=%d", names(sizes), sizes), collapse = " x "),
                 N), call. = FALSE)
  invisible(NULL)
}
