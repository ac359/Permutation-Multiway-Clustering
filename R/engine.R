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
#' for `n_reps` independent random permutation groups, aggregates the per-rep
#' p-values by `ci_agg` (\code{\link{.agg_pvals}}), and (optionally) inverts
#' the test to a confidence set. All six front ends differ only in how they
#' build their permutations and validate their inputs; everything downstream of
#' that is handled here.
#'
#' The aggregation rule is applied in exactly one place for the p-value and one
#' for the confidence set, and it is the same rule, so the reported decision and
#' the reported set are guaranteed consistent:
#'
#'   confidence set = { b : agg_r p_r(b) > alpha }
#'   p-value        =   agg_r p_r(beta_null)
#'
#' What is REPORTED for d = 1 is the closure of that set (.exact_ci_set()), so
#' an end point may itself be rejected while every interior point is accepted;
#' the direction is outward, and nothing the test accepts is ever left out.
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
#' @param ci_agg cross-rep aggregation rule, "median" (default) or "median2";
#'   see .agg_pvals(). Applied to the reported p-value and to every
#'   confidence-set path, and it also sets the rejection floor: "median2"
#'   reports min(1, 2 * median), so its smallest attainable p-value is 2/(K+1)
#'   rather than 1/(K+1), and the coarse-resolution note and the
#'   confidence-set gate below use that doubled floor.
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
  unidentified <- sum(D_resid0 * D_resid0) <= 1e-16 * sum(D * D)
  if (unidentified)
    warning(paste0("`d` has no variation after partialling out `x` (it is ",
                   "constant or collinear with the nuisance covariates): ",
                   "beta is unidentified, the test is uninformative, and ",
                   "p = 1 by construction."), call. = FALSE)
  ## Having settled the global case, a per-permutation degeneracy is a real
  ## pathology rather than the expected answer, so .ipt_prepare() stops on it
  ## instead of silently returning a noise-decided statistic (see its
  ## `degenerate` argument).
  degenerate <- if (unidentified) "zero" else "stop"

  ## Reference OLS estimate / scale (centre + step size for the CI search)
  ref <- .ols_reference(y, D, X)

  ## res_min: the p-value resolution, i.e. the spacing of the discrete grid the
  ## per-rep p-values live on, 1/(K+1). Reported as `resolution` and printed as
  ## such.
  res_min <- 1 / (K + 1L)
  ## p_floor: the smallest value the REPORTED p-value can attain, which is the
  ## grid step only under the default aggregation. Under ci_agg = "median2" the
  ## reported value is min(1, 2 * median), so the floor is 2/(K+1). The
  ## rejection gate below must consult the floor and not the step: at K = 19 and
  ## alpha = 0.05 the step is exactly alpha (so "median" can reject) while
  ## "median2" can never return below 0.10, and gating on the step there gave an
  ## unbounded interval with no note.
  agg_mult <- if (identical(ci_agg, "median2")) 2L else 1L
  p_floor  <- min(1, agg_mult * res_min)
  ## Wording and arithmetic shared by the two resolution notes below.
  floor_lab <- if (agg_mult == 2L) "2/(K+1)" else "1/(K+1)"
  ## p_floor <= alpha requires K + 1 >= agg_mult / alpha.
  need_lvl  <- ceiling(agg_mult / alpha)
  agg_why   <- if (agg_mult == 2L)
    paste0(" The floor is twice the p-value grid step 1/(K+1) = ",
           sprintf("%.3g", res_min), " because aggregate = \"median2\" ",
           "reports min(1, 2 x median).") else ""
  ## The permuted-D projections (W in the prep object) are needed for CI / joint
  ## region inversion and, in the point p-value, whenever the null is non-zero
  ## (they carry the b-statistic slope). They can be skipped only for a no-CI
  ## test of beta = 0 (e.g. a Monte-Carlo size simulation) to save work.
  want_ci     <- isTRUE(conf_int) && d == 1L && p_floor <= alpha  # interval
  want_region <- isTRUE(conf_int) && d >  1L && p_floor <= alpha  # region
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
                         cl = psock_cl, degenerate = degenerate,
                         design = type)
    list(prep = prep, pv = .ipt_eval(prep, beta0)$pvalue)
  }
  reps <- .plapply(seeds, one_rep, n_cores = if (rep_axis) n_cores else 1L)
  prep_list <- lapply(reps, `[[`, "prep")   # cached cross products, one per rep
  pv <- vapply(reps, `[[`, numeric(1), "pv")  # per-rep p-value at the null
  ## Reported p-value: the SAME aggregation the confidence set inverts, so no
  ## value the test accepts can fall outside the reported set (the set's end
  ## points are its closure -- see .exact_ci_set()). With ci_agg = "median" (the
  ## default) and any n_reps this is exactly stats::median(pv), bit for bit.
  pvalue <- .agg_pvals(matrix(pv, nrow = 1L), ci_agg)[1L]
  Kp1 <- prep_list[[1L]]$Kp1            # realised group order (K + 1)

  ## Confidence set by test inversion: an interval for a single coefficient,
  ## a joint region (grid-based) for several.
  ci <- NULL
  conf_set <- NULL
  ci_method <- NULL
  conf_region <- NULL
  conf_box <- NULL
  note <- character(0)                  # human-readable caveats appended below
  warned_res <- FALSE                   # coarse-resolution note added once
  if (isTRUE(conf_int)) {
    if (p_floor > alpha) {
      ## "Increase K" was the old advice and is not actionable on its own: K is
      ## capped by the design, so name the concrete requirement instead.
      note <- c(note, sprintf(
        paste0("No %.0f%% confidence %s: the smallest attainable p-value is ",
               "%s = %.3g, which is above alpha = %.3g, so no value ",
               "could be excluded and the set would be the whole line. A ",
               "%.0f%% set needs K + 1 >= %d -- that is, at least %d levels ",
               "in the smallest permuted dimension. The p-value reported ",
               "above is unaffected and remains exact.%s"),
        100 * conf_level, if (d == 1L) "interval" else "region",
        floor_lab, p_floor, alpha, 100 * conf_level,
        need_lvl, need_lvl, agg_why))
      warned_res <- TRUE
    } else if (d == 1L) {
      ci <- .invert_ci(prep_list, alpha = 1 - conf_level,
                       centre = ref$estimate, scale = ref$se,
                       y = y, D = D, grid = grid, agg = ci_agg)
      conf_set <- attr(ci, "conf_set")
      ci_method <- attr(ci, "ci_method")
      if (identical(ci_method, "bisection"))
        note <- c(note, paste0(
          "Confidence interval by outward bracketing and bisection, not the ",
          "exact set: the exact inversion would have to evaluate more than ",
          "the candidate budget (about 2 K^2 x n_reps points, with K = ",
          K, " and n_reps = ", n_reps, "). The end points are accurate to the ",
          "bisection tolerance and the set is assumed connected apart from ",
          "the island guard. Raise ",
          "options(mwperm.ci_exact_budget = ) to force the exact set."))
      if (isTRUE(attr(ci, "disconnected"))) {
        note <- c(note, if (identical(ci_method, "bisection")) paste0(
          "The acceptance region of the inverted test is disconnected ",
          "(per-permutation estimates that the test accepts lie outside the ",
          "interval around the OLS estimate). The reported interval was ",
          "widened to the hull of the detected components; it is ",
          "conservative.") else sprintf(paste0(
          "The confidence set is disconnected: %d separate components (see ",
          "the `conf_set` field). `conf_int` reports their hull, which is ",
          "conservative -- values inside the hull but between components are ",
          "rejected by the test."), nrow(conf_set)))
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
                            d_names = d_names, grid = grid, agg = ci_agg)
      conf_region <- reg$points
      conf_box <- reg$box
      if (length(reg$note)) note <- c(note, reg$note)
    }
  }

  if (p_floor > alpha && !warned_res) {
    note <- c(note, sprintf(
      paste0("The smallest attainable p-value is %s = %.3g, which is ",
             "above alpha = %.3g, so this test cannot reject at that level ",
             "however strong the effect. Rejecting at alpha = %.3g needs ",
             "K + 1 >= %d -- at least %d levels in the smallest permuted ",
             "dimension. The p-value itself is still exact and valid.%s"),
      floor_lab, p_floor, alpha, alpha, need_lvl, need_lvl, agg_why))
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
      conf_set    = conf_set,      # exact components of the set (d == 1):
                                   #   2-column matrix, one row per connected
                                   #   component; conf_int is its hull
      ci_method   = ci_method,     # "exact" / "grid" / "bisection", or NULL
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
      resolution  = res_min,       # p-value grid step 1/(K+1)
      p_floor     = p_floor,       # smallest attainable REPORTED p-value:
                                   #   the grid step, doubled under "median2".
                                   #   The gate for the notes above and for the
                                   #   confidence set, and what confint() reads
      note        = note,          # character vector of caveats
      call        = call           # the originating front-end call
    ),
    class = "mwperm"
  )
}

#' Aggregate per-rep p-values into the one number the confidence set inverts.
#'
#' The confidence set is defined ONCE, in one place, as
#'
#'   { b : agg_r p_r(b) > alpha }
#'
#' and every inversion path -- the exact set (.exact_ci_set), the
#' explicit-`grid` path, the bracketing fallback, and the joint region
#' (.invert_region) -- calls this function to get that number. Before 0.3.0 the three
#' single-coefficient paths disagreed: the default path took the MEDIAN OF THE
#' PER-REP END POINTS (not an inversion of anything), and the grid path took a
#' UNION over reps (systematically wider, and growing with `n_reps`). Neither
#' inverted the same function as the reported p-value.
#'
#' @param P numeric matrix, one row per candidate `b`, one column per rep.
#' @param agg "median" (Guo, Toulis and Wang 2026, Remark 1 -- the default, and
#'   the rule the reported p-value uses), "median2" (min(1, 2 * median), a
#'   genuine level-alpha p-value under arbitrary dependence across reps;
#'   Ruschendorf 1982, Vovk and Wang 2020), or "union" (the pre-0.2.0
#'   maximum-over-reps behaviour, kept only so regression tests can reproduce
#'   it).
#' @return numeric vector, one aggregated p-value per row of `P`.
#' @keywords internal
#' @noRd
.agg_pvals <- function(P, agg = "median") {
  P <- as.matrix(P)
  if (agg == "union") return(apply(P, 1L, max))
  m <- if (ncol(P) == 1L) P[, 1L] else .row_median(P)
  if (agg == "median2") pmin(1, 2 * m) else m
}

#' Row medians of one block of rows. Helper for .row_median(); see there.
#' @keywords internal
#' @noRd
.row_median_block <- function(P) {
  r <- ncol(P)
  m <- nrow(P)
  ## One radix sort of (row, value) lays every row's order statistics out
  ## contiguously: row i occupies (i-1)*r + 1 .. i*r. order() sorts NA last
  ## WITHIN each row key, not globally, so a row containing NA does not shift
  ## the blocks of later rows (verified over 4000 random NA/NaN/Inf matrices).
  v <- as.vector(P)
  sv <- v[order(rep.int(seq_len(m), r), v, method = "radix")]
  half <- (r + 1L) %/% 2L              # the same `half` median.default() uses
  if (r %% 2L == 1L)
    return(sv[seq.int(half, by = r, length.out = m)])   # selection only: exact
  lo <- sv[seq.int(half,      by = r, length.out = m)]
  hi <- sv[seq.int(half + 1L, by = r, length.out = m)]
  ## Even r: median.default() averages the two central order statistics with
  ## mean(). Call mean() once per DISTINCT (lo, hi) pair rather than once per
  ## row, so the arithmetic is byte-for-byte the arithmetic median.default()
  ## would have done. Per-rep p-values take at most K+1 distinct values, so the
  ## pair table stays small however many candidate points there are.
  uv <- unique(c(lo, hi))
  nu <- length(uv)
  pid <- (match(lo, uv) - 1L) * nu + match(hi, uv)      # exact on doubles
  up <- unique(pid)
  val <- vapply(up, function(q)
    mean(c(uv[(q - 1L) %/% nu + 1L], uv[(q - 1L) %% nu + 1L])), numeric(1))
  val[match(pid, up)]
}

#' Row-wise median, bit-identical to apply(P, 1L, stats::median).
#'
#' The exact confidence set evaluates the aggregated p-value at every breakpoint
#' and every cell between breakpoints, so `P` has O(K^2 * n_reps) rows --
#' 60,841 x 10 for a 40 x 40 dyadic fit at defaults. apply() then makes one
#' R-level median() call per row, which profiled as ~49% of the whole fit.
#'
#' Bit-identity is a HARD requirement, not a nicety: this value is compared to
#' alpha with `>`, the per-rep p-values sit exactly on the grid j/(K+1), and
#' ties at alpha are common, so a last-bit change flips an acceptance decision
#' and moves a reported interval end point. Note that the obvious spelling
#' FAILS that test: `(lo + hi)/2` is NOT bit-identical to `mean(c(lo, hi))`,
#' because mean() accumulates in LDOUBLE and then applies a second-pass
#' correction. Hand-coding the correction would be platform-dependent too
#' (LDOUBLE is 80-bit on x86, 64-bit on arm64). Hence the two devices used
#' here: order statistics by SELECTION (pure comparison, no arithmetic), and
#' mean() itself for the even case, called once per distinct pair.
#'
#' Rows are processed in blocks purely to bound peak memory, the same reason
#' .pval_matrix() chunks; each row is independent, so the block size cannot
#' change a result.
#'
#' @param P numeric matrix with at least one column.
#' @return numeric vector of length nrow(P).
#' @keywords internal
#' @noRd
.row_median <- function(P, block = 200000L) {
  if (ncol(P) == 1L) return(P[, 1L])
  m <- nrow(P)
  if (m == 0L) return(numeric(0))
  out <- if (m <= block) .row_median_block(P) else {
    o <- numeric(m)
    for (i0 in seq.int(1L, m, by = block)) {
      ii <- i0:min(i0 + block - 1L, m)
      o[ii] <- .row_median_block(P[ii, , drop = FALSE])
    }
    o
  }
  ## median.default() returns NA for any row containing NA (na.rm = FALSE).
  if (anyNA(P)) out[rowSums(is.na(P)) > 0L] <- NA_real_
  out
}

#' Per-rep p-values at a vector of candidate null values (d = 1).
#'
#' Vectorized counterpart of \code{\link{.ipt_eval}}: the statistic is affine in
#' `b`, so a whole vector of candidates costs one K x length(b) outer product
#' per rep instead of one R-level call per candidate. Chunked so a large
#' candidate set cannot blow up memory.
#'
#' @param prep_list list of per-rep prep objects (`d == 1`, `has_perm_D`).
#' @param b numeric vector of candidate null values.
#' @return numeric matrix, length(b) rows x length(prep_list) columns.
#' @keywords internal
#' @noRd
.pval_matrix <- function(prep_list, b) {
  nb <- length(b)
  P <- matrix(0, nb, length(prep_list))
  chunk <- 20000L                       # cap the K x chunk intermediates
  for (r in seq_along(prep_list)) {
    pr <- prep_list[[r]]
    u <- pr$u[1L, ]
    M <- pr$M[1L, 1L, ]
    v <- pr$v[1L, ]
    W <- pr$W[1L, 1L, ]
    K <- pr$K
    for (i0 in seq(1L, nb, by = chunk)) {
      ii <- i0:min(i0 + chunk - 1L, nb)
      bb <- b[ii]
      ## a_j(b) = |u_j - M_j b|, b_k(b) = |v_k - W_k b|, streamed one index at
      ## a time. Each element is one multiply and one subtract on scalars, so
      ## there is no reduction whose order could change; the two reductions
      ## that do exist -- min, and the count -- are a comparison and an integer
      ## sum, both exactly order-independent. So this is bit-identical to
      ## forming the two K x chunk matrices plus rep(amin, each = K), and it
      ## never allocates them (measured 1.4-1.9x on this block).
      amin <- abs(u[1L] - M[1L] * bb)
      if (K > 1L) for (j in 2L:K) amin <- pmin(amin, abs(u[j] - M[j] * bb))
      ## same comparison and tie direction as .ipt_eval(): b_k >= min_j a_j
      cnt <- integer(length(bb))
      for (k in seq_len(K)) cnt <- cnt + (abs(v[k] - W[k] * bb) >= amin)
      P[ii, r] <- (1 + cnt) / pr$Kp1
    }
  }
  P
}

#' Breakpoints of the p-value step function (d = 1).
#'
#' For a single coefficient the Procedure 1 statistics are
#'
#'   a_j(b) = |u_j - M_j b|,   b_k(b) = |v_k - W_k b|
#'
#' both piecewise linear in b. The p-value counts how many k satisfy
#' b_k(b) >= min_j a_j(b), so it can only change where two of these lines
#' cross. Dropping the absolute values, a crossing solves
#' v_k - W_k b = +/- (u_j - M_j b), giving the two root families
#'
#'   b = (u_j - v_k) / (M_j - W_k)      [the + branch]
#'   b = (u_j + v_k) / (M_j + W_k)      [the - branch]
#'
#' over all (j, k) in 1..K. Non-finite roots (a zero denominator: the two lines
#' are parallel and never cross) are discarded. Between consecutive roots the
#' p-value is exactly constant, which is what makes the confidence set
#' computable in closed form rather than by search.
#'
#' @param prep_list list of per-rep prep objects.
#' @return sorted numeric vector of distinct finite roots, pooled over reps.
#' @keywords internal
#' @noRd
.ci_breakpoints <- function(prep_list) {
  out <- vector("list", length(prep_list))
  for (r in seq_along(prep_list)) {
    pr <- prep_list[[r]]
    u <- pr$u[1L, ]
    M <- pr$M[1L, 1L, ]
    v <- pr$v[1L, ]
    W <- pr$W[1L, 1L, ]
    b <- c(outer(u, v, `-`) / outer(M, W, `-`),   # v_k - W_k b = +(u_j - M_j b)
           outer(u, v, `+`) / outer(M, W, `+`))   # v_k - W_k b = -(u_j - M_j b)
    out[[r]] <- b[is.finite(b)]
  }
  sort(unique(unlist(out, use.names = FALSE)))
}

#' Exact confidence set by test inversion (d = 1).
#'
#' Computes { b : agg_r p_r(b) > alpha } exactly, rather than approximating it.
#' Procedure 1 step 3 defines the confidence region as that set; the aggregated
#' p-value is a step function of b whose only jumps are at the roots from
#' .ci_breakpoints(), so evaluating it at every root and at one interior point
#' of every interval between consecutive roots (plus one point beyond each end)
#' determines the set completely. No bracketing assumption, no bisection
#' tolerance, and disconnected sets are returned as they are instead of being
#' replaced by their hull.
#'
#' The components are reported CLOSED. A maximal run of accepted atoms that
#' begins or ends inside an open cell is reported with the bounding breakpoint
#' as its end point, and that breakpoint is a value the run excluded -- so a
#' returned component is the topological CLOSURE of { b : agg_r p_r(b) > alpha }
#' and not the set itself. An end point may therefore be rejected by the very
#' test this inverts, while every point strictly inside the component is
#' accepted. This is the usual convention for a discrete p-value: the acceptance
#' set is a finite union of pieces whose exact end points are generally not
#' attained, so there is no attained value to report there, and closing the
#' component errs OUTWARD -- the reported set never omits an accepted value.
#' Reporting the attained side instead is a different object and would move
#' published end points; see NEWS and CLAUDE.md 5.1b.
#'
#' @param prep_list list of per-rep prep objects (`d == 1`, `has_perm_D`).
#' @param alpha test level.
#' @param agg cross-rep aggregation rule; see .agg_pvals().
#' @param budget maximum number of candidate points to evaluate. The root count
#'   grows as 2 * K^2 * n_reps, so a large group on many reps is capped; the
#'   caller then falls back to bracketing.
#' @return `NULL` when the budget is exceeded (the caller must fall back);
#'   otherwise a two-column matrix of interval end points, one row per
#'   connected component, ordered and disjoint, each component CLOSED (see
#'   Details: an end point need not itself be accepted). A zero-row matrix
#'   means the set is empty.
#' @keywords internal
#' @noRd
.exact_ci_set <- function(prep_list, alpha, agg = "median", budget = 2e5) {
  K <- prep_list[[1L]]$K
  ## Root count before forming any of them: 2 sign branches x K x K per rep.
  if (2 * as.double(K)^2 * length(prep_list) > budget) return(NULL)

  roots <- .ci_breakpoints(prep_list)
  m <- length(roots)
  if (m == 0L) {
    ## No crossing anywhere: the p-value is constant on the whole line.
    p0 <- .agg_pvals(.pval_matrix(prep_list, 0), agg)
    return(if (p0 > alpha) matrix(c(-Inf, Inf), 1L, 2L)
           else matrix(numeric(0), 0L, 2L))
  }

  ## Candidates, interleaved so that atom 1 is the unbounded cell below the
  ## smallest root, even atoms are the roots themselves, odd atoms in between
  ## are interior points of the open cells, and the last atom is the unbounded
  ## cell above the largest root.
  mids <- if (m >= 2L) (roots[-m] + roots[-1L]) / 2 else numeric(0)
  cand <- numeric(2L * m + 1L)
  cand[1L] <- roots[1L] - (1 + abs(roots[1L]))       # strictly below all roots
  cand[2L * seq_len(m)] <- roots
  if (m >= 2L) cand[2L * seq_len(m - 1L) + 1L] <- mids
  cand[2L * m + 1L] <- roots[m] + (1 + abs(roots[m]))  # strictly above all

  ## Bounds of the region each atom stands for. A root atom is the degenerate
  ## interval [r, r]; a cell atom is the open interval between its neighbours.
  lower <- c(-Inf, rep(roots, each = 2L))
  upper <- c(rep(roots, each = 2L), Inf)

  acc <- .agg_pvals(.pval_matrix(prep_list, cand), agg) > alpha
  acc[is.na(acc)] <- FALSE             # a degenerate rep cannot leak in
  idx <- which(acc)
  if (!length(idx)) return(matrix(numeric(0), 0L, 2L))

  ## Maximal runs of adjacent accepted atoms are the connected components; an
  ## accepted breakpoint next to an accepted cell is absorbed into it.
  cuts <- c(0L, which(diff(idx) > 1L), length(idx))
  n_comp <- length(cuts) - 1L
  comp <- matrix(0, n_comp, 2L)
  for (i in seq_len(n_comp)) {
    a <- idx[cuts[i] + 1L]
    z <- idx[cuts[i + 1L]]
    comp[i, ] <- c(lower[a], upper[z])
  }
  comp
}

#' Test-inversion confidence set for a single coefficient.
#'
#' Procedure 1 step 3 of Guo, Toulis and Wang (2026) defines the confidence
#' region as { b : pval(b) > alpha }. This function computes that set.
#' Permutations are held fixed across candidate values of b, so the p-value is a
#' deterministic step function of b within each rep, and every evaluation is
#' read off the cached prep objects (.ipt_prepare) in O(K) -- the whole search
#' costs no extra matrix factorizations.
#'
#' Three paths produce the set, all inverting the SAME aggregated p-value
#' agg_r p_r(b) (see .agg_pvals), and all reporting the connected components in
#' the `"conf_set"` attribute with their hull as the returned interval:
#' \describe{
#'   \item{exact (default)}{.exact_ci_set() evaluates the step function at
#'     every breakpoint and every cell between breakpoints, giving the set with
#'     no tolerance and no connectedness assumption. Components are reported
#'     closed, so an end point can be a rejected breakpoint bounding an
#'     accepted open cell; the interior is always accepted (see
#'     .exact_ci_set()).}
#'   \item{explicit `grid`}{the retained grid points, hulled; end points are
#'     accurate to the grid spacing and a set reaching a grid edge is reported
#'     as infinite there.}
#'   \item{bracketing fallback}{used only when the exact path's candidate count
#'     would exceed its budget (roughly 2 * K^2 * n_reps): outward bracketing
#'     then bisection to `tol_factor * step`, on the aggregated p-value.
#'     Accepted per-permutation estimates outside the bracket flag a
#'     disconnected set and widen the interval to the hull.}
#' }
#'
#' @param prep_list list of per-rep prep objects (each with `has_perm_D =
#'   TRUE`).
#' @param alpha,centre,scale test level, and the OLS estimate / naive SE used
#'   only to place and scale the bracketing fallback.
#' @param y,D the (unshifted) outcome and single covariate, used only to derive
#'   a sensible step size when the naive SE is unavailable.
#' @param grid optional explicit numeric grid of candidate `b`; when supplied
#'   the interval is the hull of the retained grid points (see Details).
#' @param agg cross-rep aggregation of the p-value; see .agg_pvals().
#' @param exact_budget candidate-count cap for the exact path; 0 forces the
#'   bracketing fallback. Overridable with
#'   \code{options(mwperm.ci_exact_budget = )}.
#' @return numeric length-2 interval (the hull of the set), carrying attributes
#'   `"conf_set"` (components), `"ci_method"`, `"disconnected"`, and, in grid
#'   mode, `"truncated"` / `"grid_step"` / `"grid_limit"`. On the exact path the
#'   components, and hence the interval, are the CLOSURE of the acceptance set;
#'   the grid and bracketing paths instead report attained, accepted points
#'   (accurate to the grid spacing / bisection tolerance).
#' @keywords internal
#' @noRd
.invert_ci <- function(prep_list, alpha, centre, scale, y, D, grid = NULL,
                       agg = c("median", "median2", "union"),
                       max_expand = 60L, tol_factor = 1e-3,
                       exact_budget = getOption("mwperm.ci_exact_budget",
                                                2e5)) {
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

  ## Aggregated p-value at one candidate: the single definition of the set.
  pval_at <- function(b) .agg_pvals(.pval_matrix(prep_list, b), agg)[1L]

  ## Per-permutation FWL point estimates u_j / M_j cached in the prep object.
  ## Each has p-value 1 IN ITS OWN REP (its identity statistic a_j is exactly
  ## zero there), so it is a natural candidate to probe -- used to rescue a
  ## rejected centre and to detect disconnected acceptance regions. Under
  ## cross-rep aggregation it is a candidate, not a certificate, so its
  ## aggregated p-value is always checked before it is used.
  bhat_all <- sort(unique(unlist(lapply(prep_list, function(prep) {
    bh <- prep$u[1L, ] / prep$M[1L, 1L, ]
    bh[is.finite(bh)]
  }), use.names = FALSE)))

  ## ---- explicit grid --------------------------------------------------------
  if (!is.null(grid)) {
    ## A candidate b is retained when the aggregated p-value exceeds alpha --
    ## the same rule as the reported p-value, the exact path and
    ## .invert_region(). The interval is the hull of the retained grid points;
    ## the components are reported in "conf_set", holes are flagged, and a
    ## retained set touching a grid edge yields an infinite limit on that side
    ## rather than a silently truncated finite one.
    g <- sort(unique(as.numeric(grid)))
    g <- g[is.finite(g)]
    if (length(g) < 2L)
      stop("`grid` must contain at least two distinct finite values.",
           call. = FALSE)

    p_agg <- .agg_pvals(.pval_matrix(prep_list, g), agg)
    acc <- !is.na(p_agg) & p_agg > alpha
    if (!any(acc)) {
      ci <- c(NA_real_, NA_real_)
      attr(ci, "conf_set")  <- matrix(numeric(0), 0L, 2L)
      attr(ci, "ci_method") <- "grid"
      return(ci)
    }

    idx  <- which(acc)
    i_lo <- idx[1L]
    i_hi <- idx[length(idx)]
    ## An acceptance set that reaches a grid edge is unbounded on that side as
    ## far as the grid can tell: report Inf there, but keep the finite outermost
    ## retained point in "grid_limit" so the caller can recover it if wanted.
    ci <- c(if (i_lo == 1L)        -Inf else g[i_lo],
            if (i_hi == length(g))  Inf else g[i_hi])
    ## components of the retained set, on the grid's own resolution
    cuts <- c(0L, which(diff(idx) > 1L), length(idx))
    cs <- matrix(0, length(cuts) - 1L, 2L)
    for (i in seq_len(nrow(cs)))
      cs[i, ] <- c(g[idx[cuts[i] + 1L]], g[idx[cuts[i + 1L]]])
    if (i_lo == 1L)          cs[1L, 1L]            <- -Inf
    if (i_hi == length(g))   cs[nrow(cs), 2L]      <-  Inf
    attr(ci, "conf_set")     <- cs
    attr(ci, "ci_method")    <- "grid"
    attr(ci, "disconnected") <- nrow(cs) > 1L
    attr(ci, "truncated")    <- c(i_lo == 1L, i_hi == length(g))
    attr(ci, "grid_step")    <- max(diff(g))
    attr(ci, "grid_limit")   <- c(g[i_lo], g[i_hi])
    return(ci)
  }

  ## ---- exact set ------------------------------------------------------------
  cs <- .exact_ci_set(prep_list, alpha = alpha, agg = agg,
                      budget = exact_budget)
  if (!is.null(cs)) {
    ci <- if (nrow(cs) == 0L) c(NA_real_, NA_real_)
          else c(min(cs[, 1L]), max(cs[, 2L]))   # hull, for compatibility
    attr(ci, "conf_set")     <- cs
    attr(ci, "ci_method")    <- "exact"
    attr(ci, "disconnected") <- nrow(cs) > 1L
    return(ci)
  }

  ## ---- bracketing fallback --------------------------------------------------
  ## Locate one end point of the acceptance component containing `start`.
  ## `direction` is -1 (lower) or +1 (upper).
  one_side <- function(direction, start = centre) {
    if (pval_at(start) <= alpha) {
      ## start rejected: restart from the nearest per-permutation estimate.
      if (length(bhat_all)) {
        cand <- bhat_all[which.min(abs(bhat_all - start))]
        if (pval_at(cand) > alpha) start <- cand else return(start)
      } else return(start)
    }
    lo <- start                        # last value known to be accepted
    hi <- NA_real_                     # first value known to be rejected
    h <- step                          # current step out from the start
    for (i in seq_len(max_expand)) {   # phase 1: expand out to bracket
      cand <- start + direction * h
      if (pval_at(cand) <= alpha) {
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
      if (pval_at(mid) > alpha) lo <- mid else hi <- mid
    }
    lo                                 # accepted side of the bracket
  }

  ## Island guard. The acceptance set is usually one interval, but nothing
  ## guarantees it: p(b) counts how many permuted statistics b_k(b) dominate
  ## min_j a_j(b), and both sides are piecewise linear in b, so that count can
  ## dip below the threshold and come back. Bracketing alone would return the
  ## component containing the centre and silently drop the rest. Any accepted
  ## per-permutation estimate outside the bracketed interval certifies a
  ## disconnected set; the interval is then extended to the boundary of the
  ## outlying component (a conservative hull, never narrower) and flagged.
  lo_r <- one_side(-1)
  up_r <- one_side(+1)
  disconnected <- FALSE
  out_lo <- bhat_all[bhat_all < lo_r]
  out_hi <- bhat_all[bhat_all > up_r]
  out_lo <- out_lo[vapply(out_lo, function(b) pval_at(b) > alpha, logical(1))]
  out_hi <- out_hi[vapply(out_hi, function(b) pval_at(b) > alpha, logical(1))]
  if (length(out_lo)) {
    disconnected <- TRUE
    lo_r <- one_side(-1, start = min(out_lo))
  }
  if (length(out_hi)) {
    disconnected <- TRUE
    up_r <- one_side(+1, start = max(out_hi))
  }
  ci <- c(lo_r, up_r)
  attr(ci, "conf_set")     <- matrix(c(lo_r, up_r), 1L, 2L)
  attr(ci, "ci_method")    <- "bisection"
  attr(ci, "disconnected") <- disconnected
  ci
}

#' Joint confidence region for several coefficients by test inversion.
#'
#' Inverts the exact joint test H0: beta = b over a grid of candidate vectors
#' `b`: each retained point is one at which the AGGREGATED test (.agg_pvals --
#' the median across reps by default, the same rule the reported p-value and the
#' interval paths use) does not reject at level `alpha`, so the retained set is
#' a finite-sample valid (1 - alpha) confidence region (discretised by the
#' grid). Cheap
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
#' @param agg cross-rep aggregation of the p-value; see .agg_pvals(). The same
#'   rule the single-coefficient paths invert, so the region and the interval
#'   are the same set definition in different dimensions.
#' @return list(points, box, note): `points` is the matrix of accepted `beta`
#'   vectors, `box` a 2 x d matrix of marginal (lower, upper) extents.
#' @keywords internal
#' @noRd
.invert_region <- function(prep_list, alpha, centre, scale, d_names,
                           grid = NULL, agg = "median",
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
    ## Same acceptance rule as the single-coefficient paths: the aggregated
    ## p-value across reps (.agg_pvals), evaluated at each candidate vector.
    P <- matrix(vapply(prep_list, function(pp)
                  vapply(seq_len(nrow(G)), function(i)
                    .ipt_eval(pp, G[i, ])$pvalue, numeric(1)),
                numeric(nrow(G))), nrow = nrow(G))
    acc <- .agg_pvals(P, agg) > alpha             # retained-flag per candidate
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
#' Rep seeds are `seed + r - 1` (see .ipt_engine), and within a rep dimension
#' (or cell, or block) `j` uses `rep_seed * stride + j`. Two (rep, j) pairs
#' collide exactly when the j-range reaches `stride`: with the historical
#' stride of 1000 that meant layouts with >= 1000 occupied cells, or
#' missing/irregular designs with >= 250 blocks (which use offsets up to 4q).
#' A collision does not break validity -- each rep's test is still exact,
#' because the seed only picks the random relabelling -- but two reps then
#' share a relabelling for that cell, so they are not independent and the
#' cross-rep aggregation is averaging fewer effective draws than it thinks.
#'
#' The fix is deliberately value-preserving in two ways. The arithmetic is done
#' in double, so a large seed can never overflow to a silent NA that
#' `set.seed()` would reject with a cryptic message; and the stride is widened
#' only by callers whose j-range would actually collide, so every design that
#' was already collision-free keeps the seeds -- and hence the seeded results --
#' it always had. Widening the stride for a design that needs it DOES change
#' that design's seeded output; that is the point, and it is recorded in NEWS.
#'
#' @param rep_seed the rep-level seed, or `NULL` for the ambient RNG.
#' @param j the within-rep offset (dimension, cell, or block slot).
#' @param stride the multiplier separating consecutive rep seeds. Must exceed
#'   the largest `j` the caller will use. Leave at the default 1000 unless the
#'   design can exceed it; callers that can (layouts with many cells, block
#'   designs with many blocks) pass `max(1000L, max_j + 1L)`.
#' @keywords internal
#' @noRd
.sub_seed <- function(rep_seed, j, stride = 1000) {
  if (is.null(rep_seed)) return(NULL)
  if (any(j >= stride))
    stop(sprintf(paste0("Internal error: sub-seed offset %d reaches the ",
                        "stride %d, which would make two repetitions share a ",
                        "permutation relabelling. Please report this."),
                 max(j), stride), call. = FALSE)   # nocov
  s <- rep_seed * stride + j
  if (abs(s) > .Machine$integer.max)
    stop(sprintf(paste0("`seed` is too large for this design's rep/sub-seed ",
                        "scheme: rep seed x %d + offset must stay inside R's ",
                        "integer seed range, so use |seed| below about %.0f."),
                 stride, floor(.Machine$integer.max / stride) - 1),
         call. = FALSE)
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
