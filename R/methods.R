## S3 methods for objects of class "mwperm".

#' @keywords internal
#' @noRd
.fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-3) "< 0.001" else formatC(p, digits = 3, format = "f")
}

#' Print a multi-way permutation test
#'
#' @param x An object of class \code{"mwperm"}.
#' @param digits Number of significant digits for the estimate and interval.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.mwperm <- function(x, digits = 4L, ...) {
  cat("\nInvariant permutation test (mwperm)\n")
  cat(strrep("-", 36), "\n", sep = "")
  cat("Design       :", x$type, "\n")
  nc <- x$n_clusters
  cat("Clusters     :",
      paste(sprintf("%s=%d", names(nc), nc), collapse = ", "),
      sprintf("(%d observations)\n", x$n_obs))
  cat("Permutations :", sprintf("K = %d  (group order %d, %d rep%s)\n",
                                x$K, x$n_perm, x$n_reps,
                                if (x$n_reps == 1L) "" else "s"))
  cat("\n")

  est <- x$estimate
  has_box <- !is.null(x$conf_box)      # TRUE when a joint region (d > 1) was computed
  ## One line per coefficient: estimate, plus a CI (coef 1, scalar case) or the
  ## marginal region bracket (joint case).
  for (k in seq_along(est)) {
    nm <- x$d_names[k]
    line <- sprintf("  %-12s estimate = %s", nm,
                    formatC(est[k], digits = digits, format = "g"))
    if (k == 1L && !is.null(x$conf_int) && length(x$conf_int) == 2L) {
      line <- paste0(line, sprintf("   %.0f%% CI [%s, %s]",
                                   100 * x$conf_level,
                                   formatC(x$conf_int[1], digits = digits, format = "g"),
                                   formatC(x$conf_int[2], digits = digits, format = "g")))
    } else if (has_box) {
      line <- paste0(line, sprintf("   %.0f%% region [%s, %s]",
                                   100 * x$conf_level,
                                   formatC(x$conf_box[1, k], digits = digits, format = "g"),
                                   formatC(x$conf_box[2, k], digits = digits, format = "g")))
    }
    cat(line, "\n")
  }
  if (has_box)
    cat(sprintf("  (joint %.0f%% confidence region: %d grid points retained; ",
                100 * x$conf_level, nrow(x$conf_region)),
        "brackets show its marginal extent)\n", sep = "")

  cat("\n")
  cat(sprintf("H0: beta = %s    p-value = %s\n",
              formatC(x$beta_null, format = "g"), .fmt_p(x$pvalue)))
  decision <- if (is.na(x$pvalue)) "undetermined" else
    if (x$pvalue <= x$alpha) sprintf("reject at alpha = %.2g", x$alpha) else
      sprintf("do not reject at alpha = %.2g", x$alpha)
  cat("Decision     :", decision, "\n")

  if (length(x$note)) {
    cat("\nNotes:\n")
    for (n in x$note) cat(strwrap(n, prefix = "  - ", width = 0.9 * getOption("width", 80)), sep = "\n")
    cat("\n")
  } else cat("\n")
  invisible(x)
}

#' Summarise a multi-way permutation test
#'
#' Returns (invisibly, after printing) a one-row-per-coefficient data frame with
#' the point estimate, naive standard error, the confidence limits (the inverted
#' interval for a single coefficient, or the marginal extent of the joint
#' confidence region for several), and the permutation p-value for
#' \eqn{H_0:\beta = b}.
#'
#' @param object An object of class \code{"mwperm"}.
#' @param ... Ignored.
#' @return A data frame, invisibly.
#' @export
summary.mwperm <- function(object, ...) {
  d <- length(object$estimate)         # number of coefficients
  ## Confidence limits per coefficient: a scalar interval populates only row 1,
  ## a joint region contributes the marginal box extents for every coefficient.
  ci_lo <- rep(NA_real_, d); ci_hi <- rep(NA_real_, d)
  if (!is.null(object$conf_int) && length(object$conf_int) == 2L) {
    ci_lo[1] <- object$conf_int[1]; ci_hi[1] <- object$conf_int[2]
  } else if (!is.null(object$conf_box)) {            # joint region: marginal extents
    ci_lo <- object$conf_box[1, ]; ci_hi <- object$conf_box[2, ]
  }
  tab <- data.frame(
    term       = object$d_names,
    estimate   = as.numeric(object$estimate),
    se_naive   = as.numeric(object$se_naive),
    conf_low   = ci_lo,
    conf_high  = ci_hi,
    p_value    = object$pvalue,
    row.names  = NULL,
    stringsAsFactors = FALSE
  )
  print.mwperm(object)
  invisible(tab)
}

#' Confidence interval from an inverted permutation test
#'
#' Extracts the test-inversion confidence set stored in a \code{"mwperm"} object:
#' the interval for a single coefficient, or the marginal extent of the joint
#' confidence region for several. Note that \code{level} is fixed at fitting time
#' (\code{1 - alpha}); a different level requires refitting with the
#' corresponding \code{alpha}.
#'
#' @param object An object of class \code{"mwperm"}.
#' @param parm Ignored; kept for generic compatibility.
#' @param level Confidence level; must match the level used at fitting, otherwise
#'   an error is raised (the interval cannot be re-derived without the stored
#'   permutations).
#' @param ... Ignored.
#' @return A matrix with the lower and upper limits, one row per coefficient. For
#'   a single coefficient this is the inverted-test interval; for several it is
#'   the \emph{marginal} extent of the joint confidence region (see
#'   \code{object$conf_region} for the full set of retained vectors).
#' @export
confint.mwperm <- function(object, parm, level = NULL, ...) {
  has_int <- !is.null(object$conf_int)
  has_box <- !is.null(object$conf_box)
  if (!has_int && !has_box)
    stop("No confidence set is stored in this object (it was not requested or the ",
         "resolution was too coarse). Refit with conf_int = TRUE.", call. = FALSE)
  if (!is.null(level) && !isTRUE(all.equal(level, object$conf_level)))
    stop(sprintf(paste0("This object stores a %.0f%% set; `level = %g` would need ",
                        "refitting with alpha = %g (the permutations are fixed at fit time)."),
                 100 * object$conf_level, level, 1 - level), call. = FALSE)
  ## Two-sided percentile column labels, e.g. "2.5 %" / "97.5 %" for a 95% set.
  pct <- 100 * c((1 - object$conf_level) / 2, 1 - (1 - object$conf_level) / 2)
  lab <- paste(format(pct, trim = TRUE), "%")
  if (has_int) {
    matrix(object$conf_int, nrow = 1, dimnames = list(object$d_names[1], lab))
  } else {
    out <- t(object$conf_box)                       # d x 2 (lower, upper)
    dimnames(out) <- list(object$d_names, lab)
    out
  }
}

#' Visualise the permutation reference distribution
#'
#' Plots the per-replication p-values when several replications were run, with
#' the reported (median) p-value and the level \code{alpha} marked; for a single
#' replication it draws a simple bar of the p-value against \code{alpha}. This
#' gives a quick read on the Monte-Carlo stability of the randomised test.
#'
#' @param x An object of class \code{"mwperm"}.
#' @param ... Passed to the underlying plotting calls.
#' @return \code{x}, invisibly.
#' @export
plot.mwperm <- function(x, ...) {
  pv <- x$pvalues_rep
  if (length(pv) >= 5L) {
    graphics::hist(pv, breaks = "Sturges",
                   xlim = c(0, max(1, max(pv))),
                   main = "mwperm: p-value across replications",
                   xlab = "p-value", col = "grey85", border = "white", ...)
    graphics::abline(v = x$pvalue, col = "steelblue", lwd = 2)
    graphics::abline(v = x$alpha, col = "firebrick", lwd = 2, lty = 2)
    graphics::legend("topright", bty = "n",
                     legend = c(sprintf("median p = %.3f", x$pvalue),
                                sprintf("alpha = %.2g", x$alpha)),
                     col = c("steelblue", "firebrick"), lwd = 2, lty = c(1, 2))
  } else {
    graphics::barplot(x$pvalue, ylim = c(0, max(1, x$pvalue, x$alpha * 4)),
                      names.arg = sprintf("H0: beta = %g", x$beta_null),
                      ylab = "p-value", col = "grey85", border = "white",
                      main = "mwperm: permutation p-value", ...)
    graphics::abline(h = x$alpha, col = "firebrick", lwd = 2, lty = 2)
    graphics::legend("topright", bty = "n", legend = sprintf("alpha = %.2g", x$alpha),
                     col = "firebrick", lwd = 2, lty = 2)
  }
  invisible(x)
}
