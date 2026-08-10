## Shared internals for all mwperm_* front ends. Not exported.

#' Parallel lapply with a serial default and a Windows PSOCK fallback.
#'
#' Runs \code{lapply(X, FUN)} on \code{n_cores} workers: forked
#' \code{parallel::mclapply} on Unix, a PSOCK cluster elsewhere (or when
#' forced via \code{method}, used by the tests). A pre-made cluster can be
#' supplied via \code{cl} and is then reused, NOT stopped -- the engine
#' creates one PSOCK cluster per fit and shares it across the per-rep calls
#' (spawning a fresh cluster per rep made parallel runs slower
#' than serial on the Windows path). With \code{n_cores = 1} it is
#' exactly \code{lapply}, so the default path stays base-R single-threaded.
#' Because every task in this package is either explicitly seeded or free of
#' RNG use, scheduling cannot perturb results: parallel output is identical to
#' serial (asserted by tests). Worker errors are re-thrown in the parent. If a
#' multithreaded BLAS is in use, R-level parallelism can oversubscribe cores;
#' where RhpcBLASctl is installed, workers pin BLAS to one thread.
#' \code{n_cores} beyond the detected core count is clamped silently here
#' (the engine warns once per fit, naming the argument).
#' @keywords internal
#' @noRd
.plapply <- function(X, FUN, n_cores = 1L, method = c("auto", "fork", "psock"),
                     cl = NULL) {
  method <- match.arg(method)
  n_cores <- suppressWarnings(as.integer(n_cores))
  if (length(n_cores) != 1L || is.na(n_cores))
    stop("`n_cores` must be a single integer >= 1.", call. = FALSE)
  n_cores <- max(1L, min(n_cores, .n_cores_max()))
  if (is.null(cl) && (n_cores == 1L || length(X) < 2L)) return(lapply(X, FUN))
  wrap <- function(x) {                # run in the worker
    if (requireNamespace("RhpcBLASctl", quietly = TRUE))
      try(RhpcBLASctl::blas_set_num_threads(1L), silent = TRUE)
    FUN(x)
  }
  use_fork <- (method == "fork") ||
    (method == "auto" && .Platform$OS.type == "unix")
  out <- if (!is.null(cl)) {
    parallel::parLapply(cl, X, wrap)   # caller owns the cluster's lifetime
  } else if (use_fork) {
    parallel::mclapply(X, wrap, mc.cores = n_cores)
  } else {
    cl1 <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(cl1), add = TRUE)
    parallel::parLapply(cl1, X, wrap)
  }
  bad <- vapply(out, inherits, logical(1), what = "try-error")
  if (any(bad)) stop(attr(out[[which(bad)[1L]]], "condition"))
  ## mclapply returns NULL (with only a warning) for tasks whose forked worker
  ## died without reporting an error (e.g. killed by the OS); fail loudly here
  ## rather than letting the NULL surface as an unrelated error downstream.
  if (any(vapply(out, is.null, logical(1))))
    stop("A parallel worker died without returning a result; ",
         "rerun with n_cores = 1 to see the underlying error.", call. = FALSE)
  out
}

#' Upper bound on worker count: the detected core count, further capped by
#' \code{getOption("mc.cores")} when the user (or R CMD check) has set it, so
#' the package honours the standard throttle. Inf when detection fails and no
#' option is set, so the clamp becomes a no-op rather than blocking a
#' legitimate request.
#' @keywords internal
#' @noRd
.n_cores_max <- function() {
  nc <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA)
  nc <- if (is.na(nc) || nc < 1L) Inf else nc
  opt <- suppressWarnings(as.integer(getOption("mc.cores")))
  if (length(opt) == 1L && !is.na(opt) && opt >= 1L) min(nc, opt) else nc
}

#' OLS reference estimate of the coefficient(s) of interest and a naive SE
#' (used only as a centre/scale for the confidence-interval search).
#' @keywords internal
#' @noRd
.ols_reference <- function(y, D, X) {
  D <- as.matrix(D)
  X <- as.matrix(X)
  d <- ncol(D)                         # number of coefficients of interest
  p <- ncol(X)                         # nuisance columns (incl. intercept)
  W <- cbind(X, D)                     # full design [nuisance | interest]
  fit <- stats::lm.fit(W, y)
  cf <- fit$coefficients
  idx <- (p + 1L):(p + d)              # D coefficient positions in `cf`
  est <- cf[idx]                       # OLS point estimate(s) of beta
  ## naive homoskedastic SE -- used only as a centre/scale for the CI search,
  ## never reported as an inferential quantity.
  res <- fit$residuals
  dfres <- length(y) - fit$rank        # residual degrees of freedom
  se <- rep(NA_real_, d)
  if (dfres > 0L) {
    sigma2 <- sum(res^2) / dfres       # homoskedastic error-variance estimate
    ## (W'W)^{-1}, or NULL if singular
    XtXi <- tryCatch(solve(crossprod(W)), error = function(e) NULL)
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
                        type, d_names, n_clusters, call, n_cores = 1L,
                        ci_agg = "median") {
  y <- as.numeric(y)
  D <- as.matrix(D)
  X <- as.matrix(X)
  N <- length(y)                     # number of observations
  d <- ncol(D)                       # number of coefficients of interest
  p <- ncol(X)                       # nuisance columns (incl. intercept)
  ## Validate before any linear algebra so bad input fails with its own name,
  ## not as an NA/coercion error deep inside lm.fit or a QR decomposition.
  .check_finite(list(y = y, d = D, x = X))
  if (!(is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
        alpha > 0 && alpha < 1))
    stop("`alpha` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  if (!(is.numeric(n_reps) && length(n_reps) == 1L && is.finite(n_reps) &&
        n_reps >= 1 && n_reps == trunc(n_reps)))
    stop("`n_reps` must be a single integer >= 1.", call. = FALSE)
  if (!(is.numeric(beta_null) && all(is.finite(beta_null))))
    stop("`beta_null` must be numeric and finite.", call. = FALSE)
  if (!length(beta_null) %in% c(1L, d))
    stop(sprintf(
      "`beta_null` must have length 1 or %d (one per column of `d`).", d),
      call. = FALSE)
  if (!is.null(grid) &&
      !all(vapply(if (is.list(grid)) grid else list(grid),
                  function(g) is.numeric(g) && length(g) >= 1L &&
                    all(is.finite(g)),
                  logical(1))))
    stop("`grid` must contain only finite numeric values.", call. = FALSE)
  n_cores <- suppressWarnings(as.integer(n_cores))
  if (length(n_cores) != 1L || is.na(n_cores))
    stop("`n_cores` must be a single integer >= 1.", call. = FALSE)
  nc_max <- .n_cores_max()
  if (n_cores > nc_max) {
    warning(sprintf(paste0("`n_cores` = %d exceeds the %d available cores ",
                           "(detectCores(), or getOption(\"mc.cores\") if ",
                           "set); using %d."),
                    n_cores, nc_max, nc_max), call. = FALSE)
    n_cores <- as.integer(nc_max)
  }
  if (!is.null(seed) &&
      !(is.numeric(seed) && length(seed) == 1L && is.finite(seed)))
    stop("`seed` must be NULL or a single finite number.", call. = FALSE)
  conf_level <- 1 - alpha            # always the complement of the test level
  ## Need more observations than the stacked projection [X | X_k] consumes; with
  ## p nuisance columns that projection has up to 2p columns.
  if (N <= 2L * p) {
    stop(sprintf(paste0("Need N > 2p for the projection to exist: N = %d, ",
                        "p = %d nuisance columns (incl. intercept). Drop ",
                        "covariates or add data."), N, p), call. = FALSE)
  }

  ## A `d` with no variation after partialling out X (constant, or collinear
  ## with a nuisance column) leaves beta unidentified: every residualized
  ## statistic is exactly zero in exact arithmetic, so p = 1 by the
  ## minorization (.ipt_prepare zeroes the noise-level slices to enforce
  ## that). Warn once, up front, so the p = 1 is not mistaken for evidence
  ## span(X) is contained in every span[X | X_k], so this check
  ## catches the global case; per-permutation degeneracy is handled silently
  ## by the slice floor.
  D_resid0 <- qr.resid(qr(X), D)
  if (sum(D_resid0 * D_resid0) <= 1e-16 * sum(D * D))
    warning(paste0("`d` has no variation after partialling out `x` (it is ",
                   "constant or collinear with the nuisance covariates): ",
                   "beta is unidentified, the test is uninformative, and ",
                   "p = 1 by construction."), call. = FALSE)

  ## Reference OLS estimate / scale (centre + step size for the CI search)
  ref <- .ols_reference(y, D, X)

  ## res_min: the p-value resolution, i.e. the smallest attainable p-value
  ## 1/(K+1). If it already exceeds alpha no confidence set can exclude
  ## anything.
  res_min <- 1 / (K + 1L)
  ## The permuted-D projections (W in the prep object) are needed for CI / joint
  ## region inversion and, in the point p-value, whenever the null is non-zero
  ## (they carry the b-statistic slope). They can be skipped only for a no-CI
  ## test of beta = 0 (e.g. a Monte-Carlo size simulation) to save work.
  want_ci     <- isTRUE(conf_int) && d == 1L && res_min <= alpha  # interval
  want_region <- isTRUE(conf_int) && d >  1L && res_min <= alpha  # region
  beta0 <- rep(beta_null, length.out = d)  # null value(s), recycled to d
  need_W <- want_ci || want_region || any(beta0 != 0)

  ## Build + prepare the permutations once per rep, then reuse the cached prep
  ## objects for the CI search (so the QR work is never repeated).
  ## seeds: one per rep (NULL = use the ambient RNG, no reproducibility).
  seeds <- if (is.null(seed)) rep(list(NULL), n_reps)
           else as.list(seed + seq_len(n_reps) - 1L)
  ## Parallel axis (n_cores > 1): the rep loop when it offers several
  ## explicitly seeded tasks -- each worker rebuilds its permutations from its
  ## own seed, so scheduling cannot change any draw -- otherwise the K loop
  ## inside .ipt_prepare (which uses no RNG at all). With seed = NULL the rep
  ## loop MUST stay serial: forked workers would inherit identical RNG states
  ## and silently duplicate the permutation draws across reps.
  rep_axis <- n_cores > 1L && n_reps > 1L && !is.null(seed)
  ## On non-fork platforms the K-axis path used to spawn a fresh PSOCK
  ## cluster inside every rep's .ipt_prepare (with the default
  ## n_reps that made parallel runs 3x SLOWER than serial on Windows).
  ## Create one cluster per fit and share it across the per-rep calls.
  psock_cl <- NULL
  if (!rep_axis && n_cores > 1L && K >= 2L && .Platform$OS.type != "unix") {
    psock_cl <- parallel::makePSOCKcluster(n_cores)
    on.exit(parallel::stopCluster(psock_cl), add = TRUE)
  }
  one_rep <- function(s) {
    op <- perm_builder(s)               # K+1 gather-vectors for this rep
    prep <- .ipt_prepare(y, D, X, op, need_perm_D = need_W,
                         n_cores = if (rep_axis) 1L else n_cores,
                         cl = psock_cl)
    list(prep = prep, pv = .ipt_eval(prep, beta0)$pvalue)
  }
  reps <- .plapply(seeds, one_rep, n_cores = if (rep_axis) n_cores else 1L)
  prep_list <- lapply(reps, `[[`, "prep")   # cached cross products, one per rep
  pv <- vapply(reps, `[[`, numeric(1), "pv")  # per-rep p-value at the null
  pvalue <- stats::median(pv)           # reported p-value: median across reps
  Kp1 <- prep_list[[1L]]$Kp1            # realised group order (K + 1)

  ## Confidence set by test inversion: an interval for a single coefficient,
  ## a joint region (grid-based) for several.
  ci <- NULL
  conf_region <- NULL
  conf_box <- NULL
  note <- character(0)                  # human-readable caveats appended below
  warned_res <- FALSE                   # coarse-resolution note added once
  if (isTRUE(conf_int)) {
    if (res_min > alpha) {
      ## "Increase K" was the old advice and is not actionable on its own: K is
      ## capped by the design, so name the concrete requirement instead.
      note <- c(note, sprintf(
        paste0("No %.0f%% confidence %s: the smallest attainable p-value is ",
               "1/(K+1) = %.3g, which is above alpha = %.3g, so no value ",
               "could be excluded and the set would be the whole line. A ",
               "%.0f%% set needs K + 1 >= %d -- that is, at least %d levels ",
               "in the smallest permuted dimension. The p-value reported ",
               "above is unaffected and remains exact."),
        100 * conf_level, if (d == 1L) "interval" else "region",
        res_min, alpha, 100 * conf_level,
        ceiling(1 / alpha), ceiling(1 / alpha)))
      warned_res <- TRUE
    } else if (d == 1L) {
      ci <- .invert_ci(prep_list, alpha = 1 - conf_level,
                       centre = ref$estimate, scale = ref$se,
                       y = y, D = D, grid = grid, agg = ci_agg)
      if (isTRUE(attr(ci, "disconnected"))) {
        note <- c(note, paste0(
          "The acceptance region of the inverted test is disconnected ",
          "(per-permutation estimates with p-value 1 lie outside the interval ",
          "around the OLS estimate). The reported interval was widened to the ",
          "hull of the detected components; it is conservative."))
      }
      ## Grid-mode caveats: an acceptance region that runs off an edge of the
      ## supplied `grid` is reported as infinite there (the grid cannot certify
      ## a finite bound), and the finite limits are only accurate to the grid
      ## spacing.
      trunc <- attr(ci, "truncated")
      if (!is.null(trunc) && any(trunc))
        note <- c(note, sprintf(paste0(
          "The acceptance region reaches the %s of the supplied `grid`, so ",
          "%s reported as infinite: the grid cannot certify a finite limit ",
          "beyond its own extent. Widen `grid` to bound %s."),
          paste(c("lower end", "upper end")[trunc], collapse = " and the "),
          if (sum(trunc) > 1L) "both are" else "that end is",
          if (sum(trunc) > 1L) "them" else "it"))
      gs <- attr(ci, "grid_step")
      if (!is.null(gs))
        note <- c(note, sprintf(paste0(
          "Grid-mode interval: end points are the outermost retained grid ",
          "points, so each limit is accurate to the grid spacing (%.3g) and is ",
          "conservative inward by at most that much."), gs))
      ci <- as.numeric(ci)             # drop the internal flag attributes
    } else {
      reg <- .invert_region(prep_list, alpha = 1 - conf_level,
                            centre = ref$estimate, scale = ref$se,
                            d_names = d_names, grid = grid)
      conf_region <- reg$points
      conf_box <- reg$box
      if (length(reg$note)) note <- c(note, reg$note)
    }
  }

  if (res_min > alpha && !warned_res) {
    note <- c(note, sprintf(
      paste0("The smallest attainable p-value is 1/(K+1) = %.3g, which is ",
             "above alpha = %.3g, so this test cannot reject at that level ",
             "however strong the effect. Rejecting at alpha = %.3g needs ",
             "K + 1 >= %d -- at least %d levels in the smallest permuted ",
             "dimension. The p-value itself is still exact and valid."),
      res_min, alpha, alpha, ceiling(1 / alpha), ceiling(1 / alpha)))
  }

  ## Assemble the returned "mwperm" object. Field meanings (read by the S3
  ## methods in methods.R):
  structure(
    list(
      pvalue      = pvalue,        # reported p-value (median over reps)
      pvalues_rep = pv,            # per-rep p-values (used by plot.mwperm)
      estimate    = stats::setNames(ref$estimate,
                                    d_names),  # OLS point estimate(s)
      se_naive    = stats::setNames(ref$se,
                                    d_names),        # naive homoskedastic SE(s)
      conf_int    = ci,            # inverted interval (d == 1), else NULL
      conf_region = conf_region,   # retained beta vectors (d > 1), or NULL
      conf_box    = conf_box,      # 2 x d extent of the region, or NULL
      conf_level  = conf_level,    # 1 - alpha
      alpha       = alpha,         # test level
      beta_null   = beta_null,     # null value tested
      K           = Kp1 - 1L,      # number of non-identity permutations
      n_perm      = Kp1,           # group order (K + 1)
      n_reps      = n_reps,        # independent repetitions aggregated
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
#' bracketing then bisection. The p-value at each candidate `b` is read off the
#' cached prep objects (\code{\link{.ipt_prepare}}) in O(K), so the whole
#' search costs no extra QR decompositions.
#'
#' Both modes now de-randomise across reps by the MEDIAN, matching the reported
#' p-value, the joint-region path (\code{\link{.invert_region}}), and Remark 1
#' of Guo, Toulis and Wang (2026). The default (bracketing) mode returns the
#' median of the per-rep interval end points; the explicit-`grid` mode returns
#' the hull of \eqn{\{b : \mathrm{median}_r\, p_r(b) > \alpha\}}. The two
#' coincide whenever each rep's acceptance set is a single interval; they can
#' differ only through the grid's discretisation or a disconnected acceptance
#' region (flagged, and widened to the hull). The old explicit-`grid` behaviour
#' -- retaining a point accepted in ANY rep, i.e. inverting the MAXIMUM p-value
#' across reps (a union) -- over-covered and grew with `n_reps`; it survives
#' only as \code{agg = "union"} for regression testing (see NEWS 0.2.0).
#' @param prep_list list of per-rep prep objects (each with `has_perm_D =
#'   TRUE`).
#' @param y,D the (unshifted) outcome and single covariate, used only to derive
#'   a sensible step size when the naive SE is unavailable.
#' @param grid optional explicit numeric grid of candidate `b`; when supplied
#'   the interval is the hull of the retained grid points (see Details).
#' @param agg cross-rep aggregation of the p-value: \code{"median"} (default,
#'   Remark 1), \code{"median2"} (\eqn{\min(1, 2\,\mathrm{median})}, valid under
#'   arbitrary dependence across reps; Ruschendorf 1982, Vovk & Wang 2020), or
#'   \code{"union"} (the legacy maximum-over-reps behaviour, kept only for
#'   regression tests). Applies to the explicit-`grid` mode.
#' @keywords internal
#' @noRd
.invert_ci <- function(prep_list, alpha, centre, scale, y, D, grid = NULL,
                       agg = c("median", "median2", "union"),
                       max_expand = 60L, tol_factor = 1e-3) {
  agg <- match.arg(agg)
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

  pval_at <- function(b, prep) .ipt_eval(prep, b)$pvalue  # p at null b

  ## Per-permutation FWL point estimates u_j / M_j cached in the prep object.
  ## Each has p-value 1 (its identity statistic a_j is exactly zero there), so
  ## every finite one is a certified member of the acceptance region -- used to
  ## rescue a rejected centre and to detect disconnected acceptance regions.
  bhat_of <- function(prep) {
    bh <- prep$u[1L, ] / prep$M[1L, 1L, ]
    bh[is.finite(bh)]
  }

  ## Locate one endpoint of the acceptance interval component containing
  ## `start`, for a single rep. `direction` is -1 (lower) or +1 (upper).
  one_side <- function(prep, direction, start = centre) {
    if (pval_at(start, prep) <= alpha) {
      ## start rejected: restart from the nearest per-permutation estimate,
      ## which is always accepted (p = 1). Degenerate [centre, centre] output
      ## only remains if even that fails (numerically impossible in practice).
      bh <- bhat_of(prep)
      if (length(bh)) {
        cand <- bh[which.min(abs(bh - start))]
        if (pval_at(cand, prep) > alpha) start <- cand else return(start)
      } else return(start)
    }
    lo <- start                        # last value known to be accepted
    hi <- NA_real_                     # first value known to be rejected
    h <- step                          # current step out from the start
    for (i in seq_len(max_expand)) {   # phase 1: expand out to bracket
      cand <- start + direction * h
      if (pval_at(cand, prep) <= alpha) {
        hi <- cand
        break
      }
      lo <- cand                       # cand accepted: advance the bracket
      h <- h * 1.6                     # geometric growth: few steps needed
    }
    if (is.na(hi)) return(direction * Inf)  # never rejected: unbounded
    ## phase 2: bisect, maintaining pval(lo) > alpha >= pval(hi)
    while (abs(hi - lo) > tol) {
      mid <- (lo + hi) / 2
      if (pval_at(mid, prep) > alpha) lo <- mid else hi <- mid
    }
    lo                                 # accepted side of the bracket
  }

  if (!is.null(grid)) {
    ## Explicit grid mode. A candidate b is retained when the ACROSS-REP
    ## aggregated p-value exceeds alpha, using the same de-randomisation rule
    ## as the reported p-value and as .invert_region(): the median over reps
    ## (Guo, Toulis & Wang 2026, Remark 1). The interval is the hull of the
    ## retained grid points; holes (disconnected acceptance) are flagged, and a
    ## retained set touching a grid edge yields an infinite limit on that side
    ## rather than a silently truncated finite one.
    g <- sort(unique(as.numeric(grid)))
    g <- g[is.finite(g)]
    if (length(g) < 2L)
      stop("`grid` must contain at least two distinct finite values.",
           call. = FALSE)

    ## length(g) x n_reps matrix of per-rep p-values
    P <- matrix(
      vapply(prep_list, function(pp) vapply(g, pval_at, numeric(1), prep = pp),
             numeric(length(g))),
      nrow = length(g))

    p_agg <- switch(agg,
      median  = apply(P, 1L, stats::median),
      ## 2 x median: valid under arbitrary dependence across reps
      ## (Ruschendorf 1982; Vovk & Wang 2020). Never narrower than "median".
      median2 = pmin(1, 2 * apply(P, 1L, stats::median)),
      ## legacy union == max over reps; kept only for regression testing.
      union   = apply(P, 1L, max))

    acc <- !is.na(p_agg) & p_agg > alpha   # NA-safe: a degenerate rep can't leak
    if (!any(acc)) return(c(NA_real_, NA_real_))

    idx  <- which(acc)
    i_lo <- idx[1L]
    i_hi <- idx[length(idx)]
    ## An acceptance set that reaches a grid edge is unbounded on that side as
    ## far as the grid can tell: report Inf there, but keep the finite outermost
    ## retained point in "grid_limit" so the caller can recover it if wanted.
    ci <- c(if (i_lo == 1L)        -Inf else g[i_lo],
            if (i_hi == length(g))  Inf else g[i_hi])
    attr(ci, "disconnected") <- any(diff(idx) > 1L)
    attr(ci, "truncated")    <- c(i_lo == 1L, i_hi == length(g))
    attr(ci, "grid_step")    <- max(diff(g))
    attr(ci, "grid_limit")   <- c(g[i_lo], g[i_hi])
    return(ci)
  }

  ## default mode: bracket each side in every rep, report the median endpoints.
  ## Island guard. The acceptance set {b : pval(b) > alpha} is usually one
  ## interval, but nothing guarantees it: p(b) counts how many permuted
  ## statistics b_k(b) dominate min_j a_j(b), and both sides are piecewise
  ## linear in b, so that count can dip below the threshold and come back.
  ## Bracketing alone would then return the component containing the centre
  ## and silently drop the rest. Every finite u_j / M_j has p = 1, so any of
  ## them outside the bracketed interval
  ## certifies a disconnected region. In that case extend the interval to the
  ## boundary of the outlying component(s) (a conservative hull; never narrower
  ## than before) and flag it so the engine can attach a note.
  lowers <- numeric(length(prep_list))
  uppers <- numeric(length(prep_list))
  disconnected <- FALSE
  for (r in seq_along(prep_list)) {
    prep <- prep_list[[r]]
    lo_r <- one_side(prep, -1)
    up_r <- one_side(prep, +1)
    bh <- bhat_of(prep)
    out_lo <- bh[bh < lo_r]
    out_hi <- bh[bh > up_r]
    if (length(out_lo)) {
      disconnected <- TRUE
      lo_r <- one_side(prep, -1, start = min(out_lo))
    }
    if (length(out_hi)) {
      disconnected <- TRUE
      up_r <- one_side(prep, +1, start = max(out_hi))
    }
    lowers[r] <- lo_r
    uppers[r] <- up_r
  }
  ci <- c(stats::median(lowers), stats::median(uppers))
  attr(ci, "disconnected") <- disconnected
  ci
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
.invert_region <- function(prep_list, alpha, centre, scale, d_names,
                           grid = NULL,
                           n_grid = 21L, spread = 6, max_points = 2e4L,
                           max_expand = 6L) {
  d <- length(centre)                              # number of coefficients
  ck <- ifelse(is.finite(centre), centre,
               0)       # per-coord grid centre (0 if estimate missing)
  sk <- ifelse(is.finite(scale) & scale > 0, scale,
               1)  # per-coord grid scale (1 if SE missing)

  ## Evaluate acceptance over a Cartesian grid: a point b is retained when its
  ## median p-value across reps exceeds alpha. `axes` is a list of d coordinate
  ## value-vectors; returns NULL (caller treats as "too large") if the product
  ## grid would exceed max_points.
  eval_grid <- function(axes) {
    if (prod(vapply(axes, length,
                    numeric(1))) > max_points) return(NULL)  # too big
    G <- as.matrix(expand.grid(axes))    # candidate beta vectors (rows)
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
      min(pts[, k]) <= min(res$axes[[k]]) || max(pts[,
                                                     k]) >= max(res$axes[[k]]),
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
        stop(sprintf(
          "`grid` list must have one vector per coefficient (d = %d).", d),
          call. = FALSE)
      grid
    } else rep(list(grid), d)
    res <- eval_grid(axes)
    if (is.null(res))
      return(list(points = NULL, box = NULL, note = sprintf(
        "Joint confidence region skipped: `grid` has > %d points.",
        max_points)))
  }

  if (!any(res$acc))
    return(list(points = res$G[0, , drop = FALSE], box = NULL,
                note = paste0("Joint confidence region is empty on the ",
                              "searched grid; pass a finer/shifted `grid`.")))

  pts <- res$G[res$acc, , drop = FALSE]            # retained beta vectors
  colnames(pts) <- d_names
  ## box: per-coordinate (lower, upper) extent of the retained set
  box <- rbind(apply(pts, 2L, min), apply(pts, 2L, max))
  dimnames(box) <- list(c("lower", "upper"), d_names)
  note <- if (touches(res))
    paste0("Joint confidence region still reaches the grid boundary (it may ",
           "be unbounded); pass an explicit `grid` to widen it.")
  else character(0)
  list(points = pts, box = box, note = note)
}

#' Validate and coerce a cluster id vector to dense 1-based integers.
#'
#' @param x the cluster id vector (any type coercible by \code{factor}).
#' @param what the user-facing argument name, used in the error message.
#' @keywords internal
#' @noRd
.dense_id <- function(x, what = "index") {
  f <- as.integer(factor(x))
  if (anyNA(f))
    stop(sprintf(paste0("`%s` contains missing values (NA); cluster ",
                        "identifiers must be complete."), what),
         call. = FALSE)
  f
}

#' Coerce the outcome to numeric, refusing factors.
#'
#' \code{as.numeric(factor)} yields the internal level codes -- silent data
#' corruption for an outcome. \code{d}/\code{x} are protected by
#' matrix coercion (their mode stays character and \code{.check_finite}
#' rejects it); \code{y} needs this explicit guard because factors are
#' numeric-coercible.
#' @keywords internal
#' @noRd
.check_y <- function(y) {
  if (is.factor(y))
    stop(paste0("`y` is a factor; the outcome must be numeric. Factors are ",
                "not coerced to their level codes -- if the labels are ",
                "numbers, convert explicitly with ",
                "as.numeric(as.character(y))."), call. = FALSE)
  as.numeric(y)
}

#' Error unless every supplied vector/matrix is numeric (or logical) with all
#' entries finite. NULLs are skipped; names label the user-facing arguments in
#' the error messages.
#' @keywords internal
#' @noRd
.check_finite <- function(vars) {
  for (nm in names(vars)) {
    v <- vars[[nm]]
    if (is.null(v) || length(v) == 0L) next
    if (!is.numeric(v) && !is.logical(v))
      stop(sprintf(paste0("`%s` must be numeric (got mode \"%s\"). Convert ",
                          "factors/characters to numeric columns first."),
                   nm, mode(v)), call. = FALSE)
    if (!all(is.finite(v)))
      stop(sprintf(paste0("`%s` contains missing or non-finite values ",
                          "(NA/NaN/Inf); mwperm requires complete data. Drop ",
                          "or impute the affected rows first."),
           nm), call. = FALSE)
  }
  invisible(NULL)
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
    if (nrow(X) != N) stop("`x` must have the same number of rows as `y`.",
                           call. = FALSE)
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
    if (!(is.numeric(K) && length(K) == 1L && is.finite(K) && K == trunc(K)
          && K >= 1))
      stop("`K` must be a single integer >= 1 (or NULL for the default).",
           call. = FALSE)
    K <- as.integer(K)
    if (K + 1L > smallest)
      stop(sprintf(paste0("The permutation group has order K + 1 = %d, but ",
                          "the smallest permuted dimension has only %d ",
                          "levels; the group cannot be larger than that. Use ",
                          "K <= %d, or leave K = NULL for the largest group ",
                          "this design supports."),
                   K + 1L, smallest, smallest - 1L), call. = FALSE)
  }
  if (K < 1L)
    stop("Not enough clusters to permute (need >= 2 in each dimension).",
                   call. = FALSE)
  K
}

#' Derive a per-rep, per-dimension seed from a rep-level seed (or NULL).
#'
#' Rep seeds are `seed + r - 1` (see .ipt_engine), so two (rep, j) pairs can
#' only collide when the j-range spans >= 1000 -- i.e. layouts with >= 1000
#' occupied cells or missing designs with >= 250 blocks. Even then each rep's
#' test remains valid (the seed only picks the random relabelling); only the
#' cross-rep independence of the median aggregation degrades. Deliberately
#' NOT changed: this scheme is FROZEN. Any new mixing would shift every
#' seeded result, invalidating the reference numbers the paper's tables and
#' the package's regression anchors are keyed to. Changing it requires a NEWS
#' entry, re-derived anchors and explicit sign-off.
#'
#' Computed in double so a large seed can never overflow to a silent NA;
#' values in R's integer seed range are identical to the old
#' integer arithmetic, so no seeded result changes. Beyond that range the
#' error names `seed` instead of surfacing as set.seed(NA)'s cryptic message.
#' @keywords internal
#' @noRd
.sub_seed <- function(rep_seed, j) {
  if (is.null(rep_seed)) return(NULL)
  s <- rep_seed * 1000 + j
  if (abs(s) > .Machine$integer.max)
    stop(paste0("`seed` is too large for the rep/sub-seed scheme (rep seed x ",
                "1000 + dimension offset must stay inside R's integer seed ",
                "range); use |seed| below about 2.1 million."), call. = FALSE)
  s
}

#' Column labels for the coefficient(s) of interest.
#'
#' Uses the column names of \code{D} when present, else the deparsed user
#' expression \code{fallback} (suffixed by column index when \code{D} has
#' several columns).
#'
#' \code{deparse()} yields ONE element for an ordinary symbol or short call,
#' but SEVERAL when the argument arrived as a value rather than an expression
#' -- which is exactly what \code{do.call(mwperm_dyadic, list(y, d, ...))}
#' does, a normal way to drive the package programmatically. A multi-element
#' deparse is never a usable label, and passing it to \code{setNames()} in the
#' engine errored with an opaque \code{'names' attribute [N] must be the same
#' length as the vector [1]}. Fall back to a generic \code{"d"} in that case.
#' Single-line deparses (every call that worked before) are untouched, so no
#' existing label or seeded result changes.
#' @keywords internal
#' @noRd
.coef_names <- function(D, fallback) {
  nm <- colnames(D)
  if (!is.null(nm)) return(nm)
  if (length(fallback) != 1L || is.na(fallback)) fallback <- "d"
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
    stop(sprintf("Each (%s) cell must appear at most once.", dims),
         call. = FALSE)
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
