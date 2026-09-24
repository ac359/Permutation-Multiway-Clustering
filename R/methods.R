## S3 methods for objects of class "mwperm".

#' @keywords internal
#' @noRd
.fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-3) "< 0.001" else formatC(p, digits = 3, format = "f")
}

#' Print a multi-way permutation test
#'
#' Compact display of the design, the permutation-group order and its
#' resolution, the point estimate and inverted confidence interval (when
#' available), and the permutation p-value with the reject/retain decision.
#'
#' The per-coefficient lines label their two quantities by provenance: the
#' point estimate is the *OLS* estimate (printed as `"OLS estimate"`), while
#' the confidence set comes from inverting the invariant permutation test
#' (printed as `"IPT CI"`, or `"IPT region"` for a joint confidence region).
#'
#' The `Resolution` line states that the p-value is exact but *discrete*:
#' each repetition's p-value can only take multiples of 1/(K+1), so the
#' smallest value a repetition can attain is 1/(K+1) itself, and the line then
#' names the smallest value the *reported* p-value can take -- the same
#' 1/(K+1) under the default `aggregate = "median"`, but 2/(K+1) under
#' `"median2"`, which reports `min(1, 2 * median)`. With an even number of
#' repetitions the median averages the two central values and can fall
#' between grid points; the line says so. Reading a p-value without that
#' context is the most common way to over- or under-state what the test has
#' shown -- a p equal to the floor is the strongest available evidence rather
#' than a precise number, and a design whose floor exceeds alpha cannot reject
#' at alpha however large the effect.
#'
#' With several coefficients the `H0` line marks the null as a *joint* test
#' over all of them, and the printed brackets are the marginal extent of one
#' joint confidence region -- not separate per-coefficient intervals.
#'
#' @param x An object of class `"mwperm"`.
#' @param digits Number of significant digits for the estimate and interval.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data. arXiv:2601.08610.
#' @seealso [mwperm_dyadic()], [summary.mwperm()].
#' @examples
#' data(trade_dyadic)
#' fit <- with(trade_dyadic,
#'             mwperm_dyadic(log_trade, log_dist,
#'                           x = cbind(log_gdp_i, log_gdp_j),
#'                           row = importer, col = exporter,
#'                           n_reps = 2, seed = 1))
#' print(fit)
#' @export
print.mwperm <- function(x, digits = 4L, ...) {
  ## formatC(format = "g") right-pads the non-finite and exactly-zero cases to
  ## a common width, so an unbounded limit printed as "[ -Inf,   Inf]". Trim
  ## it. Finite values are never padded, so no displayed number changes.
  fmt <- function(v) trimws(formatC(v, digits = digits, format = "g"))
  wrap <- function(txt, initial, prefix)   # notes / long explanatory lines
    cat(strwrap(txt, initial = initial, prefix = prefix,
                width = 0.9 * getOption("width", 80)), sep = "\n")

  ## A sign-flip fit relabels nothing, so its header names its own group.
  cat(if (is.null(x$n_flip)) "\nInvariant permutation test (mwperm)\n"
      else "\nInvariant sign-flip test (mwperm)\n")
  cat(strrep("-", 36), "\n", sep = "")
  cat("Design       : ", x$type, "\n", sep = "")
  if (!is.null(x$auto))                 # dispatched via mwperm(): say why
    cat("Auto-detected: ", x$auto$design, " (", x$auto$reason, ")\n", sep = "")
  nc <- x$n_clusters
  ## The layout design permutes replicates WITHIN cells, so its two counts are
  ## a cell count and the smallest cell's replicate count -- calling them
  ## "clusters" invites reading `min_cell` as a number of clusters, when it is
  ## in fact what sets K (= min_cell - 1). Every other design does permute
  ## cluster dimensions, and keeps the generic label.
  if (all(c("ncell", "min_cell") %in% names(nc))) {
    cat("Cells        : ",
        sprintf("%d cells, smallest holds %d replicates (%d observations)\n",
                nc[["ncell"]], nc[["min_cell"]], x$n_obs), sep = "")
  } else {
    cat("Clusters     : ",
        paste(sprintf("%s=%d", names(nc), nc), collapse = ", "),
        sprintf(" (%d observations)\n", x$n_obs), sep = "")
  }
  ## The sign-flip design's group is not "K permutations": its order is
  ## 2^(n_flip - 1), and n_flip -- not K -- is the number the user chose and
  ## can change. Name it that way; every permutation design keeps its line.
  if (!is.null(x$n_flip)) {
    cat("Sign flips   : ",
        sprintf("n_flip = %d flip groups  (group order 2^%d = %d, %d rep%s)\n",
                x$n_flip, x$n_flip - 1L, x$n_perm, x$n_reps,
                if (x$n_reps == 1L) "" else "s"), sep = "")
  } else {
    cat("Permutations : ", sprintf("K = %d  (group order %d, %d rep%s)\n",
                                   x$K, x$n_perm, x$n_reps,
                                   if (x$n_reps == 1L) "" else "s"), sep = "")
  }
  ## The p-value is exact but DISCRETE -- each rep's p-value can only land on
  ## multiples of 1/(K+1). Saying so here is what stops the smallest
  ## attainable value from being read as a coincidence, or a borderline one as
  ## refinable by asking for more permutations (it is not: K is capped by the
  ## design). The REPORTED value is the aggregate over reps, and the line must
  ## describe that one: its floor is `p_floor` (2/(K+1) under "median2", which
  ## reports min(1, 2 x median)), and with an even n_reps the median of the
  ## two central reps can fall between grid points. Formatted with .fmt_p, the
  ## same way the p-value below is: showing the grid step to different
  ## precision than the value sitting on it (0.0455 vs 0.045) reads as a
  ## discrepancy.
  p_floor <- if (is.null(x$p_floor)) x$resolution else x$p_floor
  cat("Resolution   : ",
      sprintf(paste0("p-values are multiples of 1/%d = %s per rep; ",
                     "reported floor %s\n"),
              x$n_perm, .fmt_p(x$resolution), .fmt_p(p_floor)), sep = "")
  if (identical(x$aggregate, "median2") ||
      (is.null(x$aggregate) && p_floor > x$resolution)) {
    wrap(paste0("aggregate = \"median2\" reports min(1, 2 x median), so the ",
                "floor is 2/(K+1)"),
         initial = "               ", prefix = "               ")
  } else if (x$n_reps > 1L && x$n_reps %% 2L == 0L) {
    wrap(sprintf("the median of %d reps can fall between grid points",
                 x$n_reps),
         initial = "               ", prefix = "               ")
  }
  cat("\n")

  est <- x$estimate
  has_box <- !is.null(x$conf_box)      # a joint region (d > 1) was computed
  ## One line per coefficient: estimate, plus a CI (coef 1, scalar case) or the
  ## marginal region bracket (joint case).
  for (k in seq_along(est)) {
    nm <- x$d_names[k]
    line <- sprintf("  %-12s OLS estimate = %s", nm, fmt(est[k]))
    if (k == 1L && !is.null(x$conf_int) && length(x$conf_int) == 2L) {
      ## An empty acceptance set is stored as [NA, NA] with a zero-row
      ## `conf_set`; print it as what it is, not as a pair of missing values.
      line <- paste0(line, if (all(is.na(x$conf_int)) &&
                                 !is.null(x$conf_set) &&
                                 nrow(x$conf_set) == 0L)
        sprintf("   %.0f%% IPT CI: empty set (see Notes)", 100 * x$conf_level)
      else sprintf("   %.0f%% IPT CI [%s, %s]", 100 * x$conf_level,
                   fmt(x$conf_int[1]), fmt(x$conf_int[2])))
    } else if (has_box) {
      line <- paste0(line, sprintf("   %.0f%% IPT region [%s, %s]",
                                   100 * x$conf_level,
                                   fmt(x$conf_box[1, k]),
                                   fmt(x$conf_box[2, k])))
    }
    cat(line, "\n", sep = "")
  }
  ## Spell out that the joint brackets are a shadow of one d-dimensional set,
  ## not d separate intervals -- reading them as per-coefficient CIs is the
  ## natural mistake, and it understates how much the coefficients trade off.
  if (has_box)
    wrap(sprintf(paste0("Joint %.0f%% confidence region over %d coefficients ",
                        "(%d grid points retained). The brackets above are ",
                        "that single region's marginal extent, not separate ",
                        "per-coefficient intervals; `conf_region` holds the ",
                        "retained vectors."),
                 100 * x$conf_level, length(est), nrow(x$conf_region)),
         initial = "  ", prefix = "  ")

  cat("\n")
  ## With several coefficients this is ONE joint null on the whole vector, and
  ## a scalar `beta_null` has been recycled to all of them -- say so, or
  ## "H0: beta = 0" reads as a single-coefficient test.
  h0 <- sprintf("H0: beta = %s",
                paste(trimws(formatC(x$beta_null, format = "g")),
                      collapse = ", "))
  if (length(est) > 1L)
    h0 <- paste0(h0, sprintf("  (joint test of %d coefficients)", length(est)))
  cat(h0, sprintf("    p-value = %s\n", .fmt_p(x$pvalue)), sep = "")
  decision <- if (is.na(x$pvalue)) "undetermined" else
    if (x$pvalue <= x$alpha) sprintf("reject at alpha = %.2g", x$alpha) else
      sprintf("do not reject at alpha = %.2g", x$alpha)
  cat("Decision     : ", decision, "\n", sep = "")

  if (length(x$note)) {
    cat("\nNotes:\n")
    for (n in x$note) wrap(n, initial = "  - ", prefix = "    ")
    cat("\n")
  } else cat("\n")
  invisible(x)
}

#' Summarise a multi-way permutation test
#'
#' Returns (invisibly, after printing) a one-row-per-coefficient data frame
#' with the OLS point estimate, naive standard error, the IPT confidence
#' limits, and the permutation p-value for H0:beta = b; see *Value* for the
#' exact meaning of each column.
#'
#' @param object An object of class `"mwperm"`.
#' @param ... Ignored.
#' @return A data frame, invisibly, with one row per coefficient and columns
#'   `term`, `ols_estimate`, `ols_se_naive`, `ipt_ci_low`, `ipt_ci_high` and
#'   `p_value`; the column names carry the provenance of each quantity.
#'   `ols_estimate` is the *OLS* point estimate (least squares on the full
#'   design) and `ols_se_naive` its naive homoskedastic OLS standard error,
#'   used internally only to centre and scale the confidence-set search -- not
#'   an inferential quantity; `ipt_ci_low`/`ipt_ci_high` are the limits of the
#'   IPT (inverted permutation test) confidence set -- the inverted-test
#'   interval for a single coefficient, or the marginal extent of the joint
#'   confidence region for several; `p_value` is the IPT permutation p-value
#'   of the joint test (repeated across rows when there are several
#'   coefficients).
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data. arXiv:2601.08610.
#' @seealso [print.mwperm()], [confint.mwperm()].
#' @examples
#' data(trade_dyadic)
#' fit <- with(trade_dyadic,
#'             mwperm_dyadic(log_trade, log_dist,
#'                           x = cbind(log_gdp_i, log_gdp_j),
#'                           row = importer, col = exporter,
#'                           n_reps = 2, seed = 1))
#' summary(fit)
#' @export
summary.mwperm <- function(object, ...) {
  d <- length(object$estimate)         # number of coefficients
  ## Confidence limits per coefficient: a scalar interval populates only row 1,
  ## a joint region contributes the marginal box extents for every coefficient.
  ci_lo <- rep(NA_real_, d)
  ci_hi <- rep(NA_real_, d)
  if (!is.null(object$conf_int) && length(object$conf_int) == 2L) {
    ci_lo[1] <- object$conf_int[1]
    ci_hi[1] <- object$conf_int[2]
  } else if (!is.null(object$conf_box)) {   # joint region: marginal extents
    ci_lo <- object$conf_box[1, ]
    ci_hi <- object$conf_box[2, ]
  }
  tab <- data.frame(
    term         = object$d_names,
    ols_estimate = as.numeric(object$estimate),
    ols_se_naive = as.numeric(object$se_naive),
    ipt_ci_low   = ci_lo,
    ipt_ci_high  = ci_hi,
    p_value      = object$pvalue,
    row.names    = NULL,
    stringsAsFactors = FALSE
  )
  print.mwperm(object)
  invisible(tab)
}

#' Confidence interval from an inverted permutation test
#'
#' Extracts the test-inversion confidence set stored in a `"mwperm"` object:
#' the interval for a single coefficient, or the marginal extent of the joint
#' confidence region for several. This is the IPT (inverted permutation test)
#' set -- obtained by inverting the finite-sample-valid test, *not* a Wald
#' interval around the OLS estimate. Note that `level` is fixed at fitting
#' time (`1 - alpha`); a different level requires refitting with the
#' corresponding `alpha`.
#'
#' @section How the set is defined: Procedure 1, step 3 of Guo, Toulis and
#'   Wang (2026) defines the confidence region as the set of null values the
#'   test does not reject, that is `{b : pval(b) > alpha}`. With `n_reps >
#'   1` there is one p-value per repetition, and the package uses a single
#'   aggregation rule everywhere. The confidence set is
#'
#' \preformatted{  { b : median over repetitions of pval(b)  >  alpha }}
#'
#' using the same median across repetitions that produces the reported p-value
#' (Remark 1 of the same paper). Under `aggregate = "median2"` the rule
#' becomes `min(1, 2 * median)` in both places. Every path computes that one
#' set: the exact inversion for a single coefficient, the explicit-`grid`
#' path, the bracketing fallback, and the joint region for several
#' coefficients. The test and the interval therefore cannot disagree in the
#' direction that matters: no value the test accepts falls outside the
#' reported set.
#'
#' For a single coefficient the set is computed exactly. The p-value is a step
#' function of the null value whose jumps are at known roots, so the package
#' evaluates it at every root and every interval between roots rather than
#' assuming the set is one interval and bisecting to a tolerance. The set need
#' not be connected; its components are in `object$conf_set` (a two-column
#' matrix of end points, one row per component), and `confint()` and
#' `object$conf_int` report their *hull*, which is conservative when there is
#' more than one component. `object$ci_method` records which path produced the
#' set: `"exact"`, `"grid"`, or `"bisection"` (the fallback used when the
#' exact candidate count would exceed its budget).
#'
#' **The reported components are the *closure* of `{b : pval(b) > alpha}`, not
#' the set itself.** The p-value is a step function, so an acceptance region
#' often begins and ends strictly between two of its jumps; there is then no
#' attained value at the boundary to report, and the exact path reports the
#' bounding jump. A printed or returned end point may therefore be a value the
#' test itself *rejects*, while every point strictly inside the interval is
#' accepted. The convention errs outward -- the interval is conservative and
#' never omits an accepted value -- and it is the reason a value exactly equal
#' to an end point should not be read as "just inside". (The `"grid"` and
#' `"bisection"` paths report attained accepted points instead, to their own
#' accuracy.)
#'
#' @param object An object of class `"mwperm"`.
#' @param parm Optional subset of coefficients: names (matching the rows of
#'   the returned matrix) or integer positions. Defaults to all coefficients.
#' @param level Confidence level; must match the level used at fitting,
#'   otherwise an error is raised (the set cannot be re-derived without the
#'   stored permutations). Defaults to the stored level.
#' @param ... Ignored.
#' @return A matrix with the lower and upper limits, one row per coefficient
#'   (rows are named by coefficient; columns keep the percentile labels the
#'   `confint` generic promises, e.g. `"2.5 \%"`/`"97.5 \%"`). The interval is
#'   the *IPT inverted-test* set, not a Wald interval: for a single
#'   coefficient the inverted-test interval, for several the *marginal* extent
#'   of the joint confidence region (see `object$conf_region` for the full set
#'   of retained vectors). The provenance is recorded in the matrix's
#'   `"method"` attribute, `"IPT (inverted permutation test)"`.
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data. arXiv:2601.08610.
#' @seealso [mwperm_dyadic()], [print.mwperm()].
#' @examples
#' data(trade_dyadic)
#' fit <- with(trade_dyadic,
#'             mwperm_dyadic(log_trade, log_dist,
#'                           x = cbind(log_gdp_i, log_gdp_j),
#'                           row = importer, col = exporter,
#'                           n_reps = 2, seed = 1))
#' confint(fit)
#' confint(fit, parm = "log_dist")
#' @export
confint.mwperm <- function(object, parm, level = NULL, ...) {
  has_int <- !is.null(object$conf_int)
  has_box <- !is.null(object$conf_box)
  if (!has_int && !has_box) {
    ## One condition, two very different remedies. When the resolution was the
    ## binding constraint, `conf_int` was already TRUE and re-passing it
    ## changes nothing -- what is needed is a design with more levels in the
    ## smallest permuted dimension. Saying "refit with conf_int = TRUE" there
    ## sends the reader down a road that cannot work.
    ## The binding quantity is the smallest attainable REPORTED p-value, which
    ## is the grid step 1/(K+1) under the default aggregation and 2/(K+1) under
    ## aggregate = "median2" (which reports min(1, 2 * median)). Objects fitted
    ## before `p_floor` existed carry only the grid step.
    floor_p <- if (is.null(object$p_floor)) object$resolution else object$p_floor
    mult <- max(1, round(floor_p * object$n_perm))   # 1, or 2 under "median2"
    ## Same vocabulary the engine's notes use: the remedy is more clusters
    ## for a permutation group, a larger n_flip for the sign-flip group.
    voc <- .group_vocab(if (is.null(object$n_flip)) "perm" else "flip",
                        as.integer(mult), ceiling(mult / object$alpha))
    why <- if (isTRUE(floor_p > object$alpha))
      sprintf(paste0("the smallest attainable p-value is %s = %.3g, ",
                     "above alpha = %.3g, so no value could have been ",
                     "excluded. That needs %s -- a refit alone will not ",
                     "produce one"),
              voc$floor, floor_p, object$alpha, voc$need)
    else "it was not requested; refit with conf_int = TRUE"
    stop(sprintf("No confidence set is stored in this object: %s.", why),
         call. = FALSE)
  }
  if (!is.null(level) && !isTRUE(all.equal(level, object$conf_level)))
    stop(sprintf(paste0("This object stores a %.0f%% set; `level = %g` ",
                        "would need refitting with alpha = %g (the ",
                        "permutations are fixed at fit time)."),
                 100 * object$conf_level, level, 1 - level), call. = FALSE)
  ## Two-sided percentile column labels, e.g. "2.5 %" / "97.5 %" for a 95% set.
  pct <- 100 * c((1 - object$conf_level) / 2, 1 - (1 - object$conf_level) / 2)
  lab <- paste(format(pct, trim = TRUE), "%")
  out <- if (has_int) {
    matrix(object$conf_int, nrow = 1, dimnames = list(object$d_names[1], lab))
  } else {
    box <- t(object$conf_box)                       # d x 2 (lower, upper)
    dimnames(box) <- list(object$d_names, lab)
    box
  }
  ## Honour the generic's subsetting contract: names or integer
  ## positions, keeping the row labels.
  if (!missing(parm) && !is.null(parm)) {
    idx <- if (is.numeric(parm)) {
      if (any(parm < 1 | parm > nrow(out) | parm != trunc(parm)))
        stop(sprintf("`parm` must index coefficients 1..%d.", nrow(out)),
             call. = FALSE)
      as.integer(parm)
    } else {
      m <- match(as.character(parm), rownames(out))
      if (anyNA(m))
        stop(sprintf("`parm` has no coefficient named %s; available: %s.",
                     paste(parm[is.na(m)], collapse = ", "),
                     paste(rownames(out), collapse = ", ")), call. = FALSE)
      m
    }
    out <- out[idx, , drop = FALSE]
  }
  ## Record provenance without touching the percentile column labels the
  ## confint generic promises (downstream code indexes by "2.5 %"/"97.5 %").
  attr(out, "method") <- "IPT (inverted permutation test)"
  out
}

## plot.mwperm lives in plot.R (with the style layer and drawing primitives).
