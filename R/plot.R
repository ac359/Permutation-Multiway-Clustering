## Plotting for "mwperm" objects. One shared style layer (.mwperm_style) feeds
## small drawing primitives (one per figure type: .mwperm_plot_coef /
## .mwperm_plot_region / .mwperm_plot_stability); plot.mwperm() dispatches on
## `type` and mwperm_save() exports at journal dimensions. Everything draws
## only from fields already stored on the object -- the test is never re-run.
## Design contract (see ?plot.mwperm): Okabe-Ito colours, but colour never
## carries information alone (marks/line types differ too, so figures survive
## grayscale printing); no top/right spines or chartjunk; every figure carries
## the design, cluster counts, N, the p-value resolution 1/(K+1), the null,
## the p-value and the decision at alpha; par() state is always restored.

#' Package plotting constants (Okabe-Ito palette, marks, lines, margins).
#'
#' Single source of truth for every figure so nothing is hard-coded at the
#' call sites. Named arguments override individual elements; unknown names
#' error (they are almost certainly typos). Grayscale/colour-vision safety:
#' the data marks (estimate, interval, region) and the reference lines (null,
#' alpha) differ in pch/lty as well as colour.
#' @keywords internal
#' @noRd
.mwperm_style <- function(...) {
  s <- list(
    ## Okabe-Ito colourblind-safe palette (Wong 2011, Nat. Methods), exposed
    ## so overrides can pick from it by name.
    palette      = c(blue = "#0072B2", vermillion = "#D55E00",
                     green = "#009E73", orange = "#E69F00",
                     skyblue = "#56B4E9", purple = "#CC79A7",
                     yellow = "#F0E442", grey = "#999999", black = "#000000"),
    ## data marks: the OLS estimate and the inverted-test interval/region
    col_estimate = "#0072B2",          # Okabe-Ito blue
    pch_estimate = 19,                 # filled circle
    cex_estimate = 1.1,
    col_interval = "#0072B2",
    lwd_interval = 2,
    cap_len      = 0.05,               # CI end-cap half-length (inches)
    col_region   = "#56B4E9",          # accepted-region points/hull
    ## reference marks: the null value and the test level alpha
    col_null     = "#D55E00",          # Okabe-Ito vermillion
    lty_null     = 2,                  # dashed: readable without colour
    lwd_null     = 1.2,
    pch_null     = 4,                  # cross, for a point null in 2-D plots
    col_alpha    = "#D55E00",
    lty_alpha    = 3,                  # dotted: never the same lty as the null
    lwd_alpha    = 1.2,
    ## fills and outlines
    col_fill     = "grey85",           # histogram bars
    col_border   = "white",
    col_box      = "grey45",           # marginal conf_box outline
    lty_box      = 3,
    ## structural ink: axes, subtitle, annotations
    col_axis     = "grey30",
    col_sub      = "grey30",
    col_annot    = "black",
    cex_main     = 1.0,
    cex_sub      = 0.72,
    cex_annot    = 0.78,
    cex_axis     = 0.85,
    cex_lab      = 0.9,
    tcl          = -0.3,               # short outward ticks
    mar          = c(4.1, 4.1, 4.1, 1.1)
  )
  ov <- list(...)
  if (length(ov)) {
    bad <- setdiff(names(ov), names(s))
    if (length(bad))
      stop(sprintf("Unknown style element(s): %s. See ?plot.mwperm for the list.",
                   paste0("`", bad, "`", collapse = ", ")), call. = FALSE)
    s[names(ov)] <- ov
  }
  s
}

#' Standard subtitle: design, per-dimension cluster counts, N, resolution
#' (plus the biclique-block summary for the missing design).
#' @keywords internal
#' @noRd
.mwperm_subtitle <- function(x) {
  nc <- x$n_clusters
  parts <- c(
    sprintf("%s design", x$type),
    paste(sprintf("%s = %d", names(nc), nc), collapse = ", "),
    if (!is.null(x$n_blocks))
      sprintf("%d block%s, %d/%d cells", x$n_blocks,
              if (x$n_blocks == 1L) "" else "s", x$cells_used, x$cells_total),
    sprintf("N = %d", x$n_obs),
    sprintf("resolution 1/(K+1) = %s",
            formatC(x$resolution, digits = 3, format = "g"))
  )
  paste(parts, collapse = " | ")
}

#' One-line statistical annotation: the null, the p-value, the decision.
#' @keywords internal
#' @noRd
.mwperm_decision <- function(x) {
  dec <- if (is.na(x$pvalue)) "undetermined" else
    if (x$pvalue <= x$alpha) sprintf("reject at alpha = %.2g", x$alpha) else
      sprintf("do not reject at alpha = %.2g", x$alpha)
  sprintf("H0: beta = %s   p = %s   (%s)",
          paste(formatC(x$beta_null, format = "g"), collapse = ", "),
          .fmt_p(x$pvalue), dec)
}

#' Shrink a cex until the text fits the panel width (never below min_cex).
#' Must be called after plot.new() (strwidth needs an open plot).
#' @keywords internal
#' @noRd
.shrink_cex <- function(txt, cex, min_cex = 0.55) {
  w <- graphics::strwidth(txt, units = "inches", cex = cex)
  avail <- graphics::par("pin")[1L]
  if (is.finite(w) && w > 0 && w > avail) max(min_cex, cex * avail / w) else cex
}

#' Shared plot skeleton: open panel, bottom/left axes only (no top/right
#' spines, outward ticks), and the three-line title block -- main, the
#' standard subtitle, and the H0 / p-value / decision annotation.
#' @keywords internal
#' @noRd
.mwperm_frame <- function(fit, style, xlim, ylim, main, sub, xlab, ylab,
                          y_axis = TRUE, yat = NULL, ylabels = NULL) {
  graphics::plot.new()
  graphics::plot.window(xlim = xlim, ylim = ylim)
  graphics::axis(1, col = style$col_axis, col.axis = style$col_axis,
                 cex.axis = style$cex_axis, tcl = style$tcl)
  if (y_axis) {
    if (is.null(yat))
      graphics::axis(2, col = style$col_axis, col.axis = style$col_axis,
                     cex.axis = style$cex_axis, tcl = style$tcl, las = 1)
    else                                # categorical rows: labels only, no
      graphics::axis(2, at = yat, labels = ylabels,  # spine/ticks joining them
                     col.axis = style$col_axis, cex.axis = style$cex_axis,
                     lwd = 0, lwd.ticks = 0, las = 1)
  }
  if (!is.null(xlab) && nzchar(xlab))
    graphics::title(xlab = xlab, line = 2.2, cex.lab = style$cex_lab,
                    col.lab = style$col_axis)
  if (!is.null(ylab) && nzchar(ylab))
    graphics::title(ylab = ylab, line = 2.9, cex.lab = style$cex_lab,
                    col.lab = style$col_axis)
  if (!is.null(main) && nzchar(main))
    graphics::mtext(main, side = 3, line = 2.5, font = 2,
                    cex = .shrink_cex(main, style$cex_main))
  if (is.null(sub)) sub <- .mwperm_subtitle(fit)
  if (nzchar(sub))
    graphics::mtext(sub, side = 3, line = 1.55, col = style$col_sub,
                    cex = .shrink_cex(sub, style$cex_sub))
  dec <- .mwperm_decision(fit)
  graphics::mtext(dec, side = 3, line = 0.55, col = style$col_annot,
                  cex = .shrink_cex(dec, style$cex_annot))
  invisible(NULL)
}

#' Pad a range for plotting limits, ignoring non-finite values.
#' @keywords internal
#' @noRd
.pad_range <- function(v, f = 0.15) {
  v <- v[is.finite(v)]
  if (!length(v)) v <- c(-1, 1)
  span <- diff(range(v))
  if (span <= 0) span <- max(abs(v), 1)
  range(v) + c(-f, f) * span
}

#' Flagship figure: OLS estimate with the inverted-test confidence set.
#'
#' d = 1: point estimate with the IPT interval (end caps; an unbounded side
#' is drawn as an outward arrow) against a dashed null reference. d > 1: a
#' forest of the joint region's marginal extents (conf_box), rows ordered by
#' estimate. Availability is checked by the dispatcher.
#' @keywords internal
#' @noRd
.mwperm_plot_coef <- function(x, style, main = NULL, sub = NULL,
                              xlab = NULL, ylab = NULL) {
  d <- length(x$estimate)
  est <- as.numeric(x$estimate)
  b0 <- rep(x$beta_null, length.out = d)

  if (d == 1L) {
    op <- graphics::par(mar = style$mar)
    on.exit(graphics::par(op), add = TRUE)
    ci <- as.numeric(x$conf_int)
    xlim <- .pad_range(c(est, b0, ci), f = 0.18)
    ylim <- c(0, 2)                    # single row at y = 1; headroom for labels
    .mwperm_frame(x, style, xlim, ylim,
                  main = if (is.null(main))
                    sprintf("OLS estimate and %.0f%% IPT confidence interval",
                            100 * x$conf_level) else main,
                  sub = sub,
                  xlab = if (is.null(xlab))
                    sprintf("coefficient on %s", x$d_names[1L]) else xlab,
                  ylab = ylab, y_axis = FALSE)
    ## null reference (dashed vertical), direct-labelled at the top
    graphics::segments(b0, 0, b0, 1.76, col = style$col_null,
                       lty = style$lty_null, lwd = style$lwd_null)
    graphics::text(b0, 1.86, "H0", col = style$col_null, cex = style$cex_annot)
    ## interval with end caps; an unbounded side becomes an outward arrow
    cap <- style$cap_len * diff(ylim) / graphics::par("pin")[2L]
    y0 <- 1
    lo <- ci[1L]; hi <- ci[2L]
    tail <- if (is.finite(est)) est else b0[1L]   # arrow tail for unbounded sides
    xl <- if (is.finite(lo)) lo else xlim[1L]
    xr <- if (is.finite(hi)) hi else xlim[2L]
    graphics::segments(xl, y0, xr, y0, col = style$col_interval,
                       lwd = style$lwd_interval, lend = 1)
    if (is.finite(lo))
      graphics::segments(lo, y0 - cap, lo, y0 + cap, col = style$col_interval,
                         lwd = style$lwd_interval, lend = 1)
    else
      graphics::arrows(tail, y0, xlim[1L], y0, length = 0.07, angle = 30,
                       col = style$col_interval, lwd = style$lwd_interval)
    if (is.finite(hi))
      graphics::segments(hi, y0 - cap, hi, y0 + cap, col = style$col_interval,
                         lwd = style$lwd_interval, lend = 1)
    else
      graphics::arrows(tail, y0, xlim[2L], y0, length = 0.07, angle = 30,
                       col = style$col_interval, lwd = style$lwd_interval)
    graphics::points(est, y0, pch = style$pch_estimate,
                     cex = style$cex_estimate, col = style$col_estimate)
    ## direct labels: interval value below the whisker, estimate above the point
    graphics::text((xl + xr) / 2, y0 - 3 * cap,
                   sprintf("%.0f%% IPT CI [%s, %s]", 100 * x$conf_level,
                           formatC(lo, digits = 3, format = "g"),
                           formatC(hi, digits = 3, format = "g")),
                   cex = style$cex_annot, col = style$col_sub)
    if (is.finite(est))
      graphics::text(est, y0 + 3 * cap,
                     sprintf("OLS estimate = %s",
                             formatC(est, digits = 3, format = "g")),
                     cex = style$cex_annot, col = style$col_sub)
  } else {
    ## forest of the joint region's marginal extents, rows ordered by estimate
    box <- x$conf_box
    ord <- order(est)                  # smallest at the bottom row
    lo <- box[1L, ord]; hi <- box[2L, ord]
    eo <- est[ord]; nn <- x$d_names[ord]; b0o <- b0[ord]
    ## widen the left margin so horizontal (las = 1) names fit
    mar <- style$mar
    mar[2L] <- max(mar[2L], min(10, 1.6 + 0.55 * max(nchar(nn))))
    op <- graphics::par(mar = mar)
    on.exit(graphics::par(op), add = TRUE)
    yy <- seq_len(d)
    xlim <- .pad_range(c(eo, b0, lo, hi), f = 0.12)
    ylim <- c(0.4, d + 0.6)
    .mwperm_frame(x, style, xlim, ylim,
                  main = if (is.null(main))
                    sprintf("OLS estimates and joint %.0f%% IPT confidence region",
                            100 * x$conf_level) else main,
                  sub = sub,
                  xlab = if (is.null(xlab)) "coefficient value" else xlab,
                  ylab = ylab, yat = yy, ylabels = nn)
    ## null reference: one dashed line for a common null, else per-row crosses
    if (length(unique(b0)) == 1L)
      graphics::segments(b0[1L], ylim[1L], b0[1L], ylim[2L],
                         col = style$col_null, lty = style$lty_null,
                         lwd = style$lwd_null)
    else
      graphics::points(b0o, yy, pch = style$pch_null, col = style$col_null,
                       lwd = 1.5, cex = 1.1)
    cap <- min(0.18, style$cap_len * diff(ylim) / graphics::par("pin")[2L])
    graphics::segments(lo, yy, hi, yy, col = style$col_interval,
                       lwd = style$lwd_interval, lend = 1)
    graphics::segments(lo, yy - cap, lo, yy + cap, col = style$col_interval,
                       lwd = style$lwd_interval, lend = 1)
    graphics::segments(hi, yy - cap, hi, yy + cap, col = style$col_interval,
                       lwd = style$lwd_interval, lend = 1)
    graphics::points(eo, yy, pch = style$pch_estimate,
                     cex = style$cex_estimate, col = style$col_estimate)
    ## honesty note: these are marginal extents of one joint region, not
    ## per-coefficient intervals
    graphics::mtext(sprintf(
      "whiskers: marginal extent of the joint %.0f%% region (not per-coefficient intervals)",
      100 * x$conf_level), side = 1, line = 3.2,
      col = style$col_sub, cex = .shrink_cex("m", style$cex_sub))
  }
  invisible(NULL)
}

#' Joint confidence region for exactly two coefficients: accepted grid points
#' with a convex-hull outline (the points stay visible, so a non-convex
#' region is not overstated by the hull), the marginal conf_box, the OLS
#' estimate and the null.
#' @keywords internal
#' @noRd
.mwperm_plot_region <- function(x, style, main = NULL, sub = NULL,
                                xlab = NULL, ylab = NULL) {
  op <- graphics::par(mar = style$mar)
  on.exit(graphics::par(op), add = TRUE)
  G <- x$conf_region
  est <- as.numeric(x$estimate)
  b0 <- rep(x$beta_null, length.out = 2L)
  xlim <- .pad_range(c(G[, 1L], est[1L], b0[1L]), f = 0.12)
  ylim <- .pad_range(c(G[, 2L], est[2L], b0[2L]), f = 0.12)
  ylim[2L] <- ylim[2L] + 0.28 * diff(ylim)   # headroom reserved for the legend
  .mwperm_frame(x, style, xlim, ylim,
                main = if (is.null(main))
                  sprintf("Joint %.0f%% IPT confidence region",
                          100 * x$conf_level) else main,
                sub = sub,
                xlab = if (is.null(xlab)) x$d_names[1L] else xlab,
                ylab = if (is.null(ylab)) x$d_names[2L] else ylab)
  ## marginal box behind everything (light, dotted)
  if (!is.null(x$conf_box))
    graphics::rect(x$conf_box[1L, 1L], x$conf_box[1L, 2L],
                   x$conf_box[2L, 1L], x$conf_box[2L, 2L],
                   border = style$col_box, lty = style$lty_box)
  ## the region: hull outline plus the accepted grid points themselves
  if (nrow(G) >= 3L) {
    h <- grDevices::chull(G[, 1L], G[, 2L])
    graphics::polygon(G[h, 1L], G[h, 2L],
                      col = grDevices::adjustcolor(style$col_region, 0.25),
                      border = style$col_region, lwd = 1.5)
  }
  graphics::points(G[, 1L], G[, 2L], pch = 15, cex = 0.4,
                   col = grDevices::adjustcolor(style$col_region, 0.8))
  graphics::points(est[1L], est[2L], pch = style$pch_estimate,
                   col = style$col_estimate, cex = style$cex_estimate)
  graphics::points(b0[1L], b0[2L], pch = style$pch_null, col = style$col_null,
                   lwd = 2, cex = 1.1)
  ## four distinct marks share the panel: a legend is unavoidable here
  graphics::legend("topright", bty = "n", cex = style$cex_annot,
                   legend = c(sprintf("accepted (%.0f%% region)", 100 * x$conf_level),
                              "marginal box", "OLS estimate", "H0 (null)"),
                   col = c(style$col_region, style$col_box,
                           style$col_estimate, style$col_null),
                   pch = c(15, NA, style$pch_estimate, style$pch_null),
                   lty = c(NA, style$lty_box, NA, NA),
                   pt.cex = c(0.7, 1, 1, 1))
  invisible(NULL)
}

#' Monte-Carlo diagnostic: stability of the randomized p-value across reps.
#'
#' With >= 5 reps, a histogram of the per-rep p-values with the reported
#' (median) p-value and alpha marked; with fewer, a one-row dot strip of the
#' per-rep p-values against alpha (degrades gracefully to a single point).
#' @keywords internal
#' @noRd
.mwperm_plot_stability <- function(x, style, main = NULL, sub = NULL,
                                   xlab = NULL, ylab = NULL) {
  op <- graphics::par(mar = style$mar)
  on.exit(graphics::par(op), add = TRUE)
  pv <- x$pvalues_rep[is.finite(x$pvalues_rep)]
  if (is.null(main))
    main <- "Diagnostic: Monte-Carlo stability of the randomized p-value"
  ## zoom to the data: an axis fixed at [0, 1] hides all detail when the
  ## p-values pile up near zero (the typical rejection case)
  xmax <- max(c(pv, x$alpha, if (is.finite(x$pvalue)) x$pvalue))
  xlim <- c(0, min(1, max(0.12, 1.35 * xmax)))
  ## direct label beside a vertical line, on the side away from the nearer
  ## panel edge
  lab <- function(at, txt, col, y)
    graphics::text(at, y, txt, col = col, cex = style$cex_annot,
                   adj = if (at > mean(xlim)) c(1.08, 0) else c(-0.08, 0))

  if (length(pv) >= 5L) {
    h <- graphics::hist(pv, breaks = "Sturges", plot = FALSE)
    ylim <- c(0, max(h$counts) * 1.2)
    .mwperm_frame(x, style, xlim, ylim, main = main, sub = sub,
                  xlab = if (is.null(xlab))
                    sprintf("per-replication p-value (n_reps = %d)", x$n_reps)
                  else xlab,
                  ylab = if (is.null(ylab)) "replications" else ylab)
    graphics::rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1L], h$counts,
                   col = style$col_fill, border = style$col_border)
    top <- ylim[2L]
    if (is.finite(x$pvalue))
      graphics::segments(x$pvalue, 0, x$pvalue, 0.9 * top,
                         col = style$col_estimate, lwd = 2)
    graphics::segments(x$alpha, 0, x$alpha, 0.9 * top, col = style$col_alpha,
                       lty = style$lty_alpha, lwd = style$lwd_alpha)
    ## labels nudged apart vertically when the two lines nearly coincide
    y_alpha <- if (is.finite(x$pvalue) &&
                   abs(x$pvalue - x$alpha) < 0.15 * diff(xlim))
      0.82 * top else 0.93 * top
    if (is.finite(x$pvalue))
      lab(x$pvalue, sprintf("median p = %s", .fmt_p(x$pvalue)),
          style$col_estimate, 0.93 * top)
    lab(x$alpha, sprintf("alpha = %.2g", x$alpha), style$col_alpha, y_alpha)
  } else {
    ## few reps: a dot strip, one point per replication
    ylim <- c(0, 2)
    .mwperm_frame(x, style, xlim, ylim, main = main, sub = sub,
                  xlab = if (is.null(xlab))
                    sprintf("p-value (%d replication%s)", x$n_reps,
                            if (x$n_reps == 1L) "" else "s")
                  else xlab,
                  ylab = ylab, y_axis = FALSE)
    graphics::segments(x$alpha, 0, x$alpha, 1.7, col = style$col_alpha,
                       lty = style$lty_alpha, lwd = style$lwd_alpha)
    lab(x$alpha, sprintf("alpha = %.2g", x$alpha), style$col_alpha, 1.75)
    if (length(pv))
      graphics::points(pv, rep(1, length(pv)), pch = style$pch_estimate,
                       col = grDevices::adjustcolor(style$col_estimate, 0.75),
                       cex = 1.1)
    if (is.finite(x$pvalue)) {
      ## keep the p label clear of the alpha line when the two nearly coincide
      near <- abs(x$pvalue - x$alpha) < 0.15 * diff(xlim)
      adj <- if (near) {
        if (x$pvalue <= x$alpha) c(1.08, 0) else c(-0.08, 0)
      } else if (x$pvalue > mean(xlim)) c(1.08, 0) else c(-0.08, 0)
      graphics::text(x$pvalue, 1.3, sprintf("p = %s", .fmt_p(x$pvalue)),
                     col = style$col_estimate, cex = style$cex_annot, adj = adj)
    }
  }
  invisible(NULL)
}

#' Plot a multi-way permutation test
#'
#' Publication-style figures for \code{"mwperm"} objects, drawn with base
#' graphics from fields stored on the object (the test is never re-run).
#' Figures use the colourblind-safe Okabe-Ito palette, and colour never
#' carries information alone -- marks and line types differ too, so the
#' figures survive grayscale printing. Every figure carries the design type,
#' per-dimension cluster counts, N, the attainable p-value resolution
#' \code{1/(K+1)}, the null value, the p-value, and the
#' reject/do-not-reject decision at \code{alpha}.
#'
#' @param x An object of class \code{"mwperm"}.
#' @param type Which figure to draw:
#'   \describe{
#'     \item{\code{"auto"} (default)}{the flagship \code{"coef"} figure whenever
#'       a confidence set is stored on the object, otherwise the
#'       \code{"stability"} diagnostic.}
#'     \item{\code{"coef"}}{the OLS point estimate against the inverted-test
#'       (IPT) confidence set: for one coefficient, the estimate with the IPT
#'       interval (end caps; an unbounded side is drawn as an outward arrow)
#'       and a dashed reference at the null; for several, a forest of the joint
#'       region's marginal extents (\code{conf_box}), rows ordered by estimate.
#'       Note the coefficient is the \emph{partialled-out} effect: nuisance
#'       covariates are projected out by the test.}
#'     \item{\code{"region"}}{the joint confidence region for exactly two
#'       coefficients: accepted grid points from \code{conf_region} with a
#'       convex-hull outline (points stay visible so a non-convex region is not
#'       overstated), the marginal box, the estimate, and the null.}
#'     \item{\code{"stability"}}{the Monte-Carlo diagnostic: per-replication
#'       p-values (histogram for \code{n_reps >= 5}, a dot strip otherwise)
#'       with the reported median p-value and \code{alpha} marked.}
#'     \item{\code{"null"}, \code{"profile"}}{reserved for the permutation null
#'       distribution and the test-inversion p-value curve; they require
#'       statistics stored at fit time, which this version does not retain, so
#'       they currently message and fall back to the default figure.}
#'     \item{\code{"all"}}{every figure available for this object, arranged on
#'       one page (\code{par("mfrow")} is restored on exit).}
#'   }
#'   A requested figure whose ingredients are not stored (e.g. \code{"coef"}
#'   after \code{conf_int = FALSE}) falls back with a message; it never errors.
#' @param ... Style overrides by name, e.g. \code{col_estimate = "black"} or
#'   \code{lwd_interval = 3}. Unknown arguments are ignored with a warning
#'   (nothing is forwarded blindly to the underlying graphics calls). The
#'   elements: \code{palette} (the Okabe-Ito colours); \code{col_estimate},
#'   \code{pch_estimate}, \code{cex_estimate}; \code{col_interval},
#'   \code{lwd_interval}, \code{cap_len}; \code{col_region}; \code{col_null},
#'   \code{lty_null}, \code{lwd_null}, \code{pch_null}; \code{col_alpha},
#'   \code{lty_alpha}, \code{lwd_alpha}; \code{col_fill}, \code{col_border};
#'   \code{col_box}, \code{lty_box}; \code{col_axis}, \code{col_sub},
#'   \code{col_annot}; \code{cex_main}, \code{cex_sub}, \code{cex_annot},
#'   \code{cex_axis}, \code{cex_lab}; \code{tcl}; \code{mar}.
#' @param main,sub,xlab,ylab Usual title overrides. \code{sub} replaces the
#'   standard subtitle line (design type, per-dimension cluster counts, N,
#'   and the p-value resolution \code{1/(K+1)}); \code{sub = ""} suppresses
#'   it.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{mwperm_save}} for export at journal dimensions;
#'   \code{\link{mwperm_dyadic}}, \code{\link{print.mwperm}}.
#' @examples
#' data(trade_dyadic)
#' fit <- with(trade_dyadic,
#'             mwperm_dyadic(y = log_trade, d = log_dist,
#'                           x = cbind(log_gdp_i, log_gdp_j),
#'                           row = importer, col = exporter,
#'                           n_reps = 7, seed = 1))
#' plot(fit)                     # flagship: OLS estimate + inverted IPT CI
#' plot(fit, type = "stability") # Monte-Carlo diagnostic across the 7 reps
#' plot(fit, type = "all")       # both, on one page
#' plot(fit, col_estimate = "black", main = "Distance elasticity")
#' @export
plot.mwperm <- function(x, type = c("auto", "coef", "region", "null",
                                    "profile", "stability", "all"), ...,
                        main = NULL, sub = NULL, xlab = NULL, ylab = NULL) {
  type <- match.arg(type)
  ## split `...`: names matching style constants become overrides; anything
  ## else is ignored with a warning rather than forwarded blindly (passing
  ## e.g. `col` used to crash hist() with "matched by multiple actual
  ## arguments").
  dots <- list(...)
  nm <- names(dots); if (is.null(nm)) nm <- rep("", length(dots))
  is_sty <- nzchar(nm) & nm %in% names(.mwperm_style())
  if (any(!is_sty) && length(dots))
    warning(sprintf("Ignoring unsupported argument(s): %s. See ?plot.mwperm for the style elements.",
                    paste0("`", nm[!is_sty], "`", collapse = ", ")),
            call. = FALSE)
  style <- do.call(.mwperm_style, dots[is_sty])

  ## what the stored fields support
  d <- length(x$estimate)
  has_ci   <- !is.null(x$conf_int) && length(x$conf_int) == 2L && !anyNA(x$conf_int)
  has_box  <- !is.null(x$conf_box)
  has_reg  <- d == 2L && !is.null(x$conf_region) && nrow(x$conf_region) >= 1L
  has_coef <- (d == 1L && has_ci) || (d > 1L && has_box)

  ## resolve "auto" and downgrade unavailable types with a message (never an
  ## error), so plot(fit) works for every design and option combination
  resolve <- function(tp) {
    if (tp %in% c("null", "profile")) {
      message(sprintf(paste0("type = \"%s\" needs permutation %s stored at fit ",
                             "time, which this object does not carry; falling ",
                             "back to the default figure."),
                      tp, if (tp == "null") "statistics" else "p-value profiles"))
      tp <- "auto"
    }
    if (tp == "auto") tp <- if (has_coef) "coef" else "stability"
    if (tp == "coef" && !has_coef) {
      message(paste0("No confidence set is stored on this object (conf_int = ",
                     "FALSE, empty inversion, or resolution 1/(K+1) > alpha); ",
                     "showing the p-value stability diagnostic instead."))
      tp <- "stability"
    }
    if (tp == "region" && !has_reg) {
      message(if (d != 2L)
        "type = \"region\" needs exactly two coefficients (d = 2); showing the default figure instead."
        else
          "No joint confidence region is stored on this object; showing the default figure instead.")
      tp <- if (has_coef) "coef" else "stability"
    }
    tp
  }

  draw <- function(tp)
    switch(tp,
           coef      = .mwperm_plot_coef(x, style, main, sub, xlab, ylab),
           region    = .mwperm_plot_region(x, style, main, sub, xlab, ylab),
           stability = .mwperm_plot_stability(x, style, main, sub, xlab, ylab))

  if (type == "all") {
    panels <- c(if (has_coef) "coef", if (has_reg) "region", "stability")
    if (length(panels) > 1L) {
      op <- graphics::par(mfrow = if (length(panels) <= 2L)
        c(1L, length(panels)) else c(2L, 2L))
      on.exit(graphics::par(op), add = TRUE)
    }
    for (tp in panels) draw(tp)
  } else {
    draw(resolve(type))
  }
  invisible(x)
}

#' Save an mwperm figure at journal dimensions
#'
#' Renders \code{\link{plot.mwperm}} to a file sized for a journal column,
#' using only base \code{grDevices} devices. The device is always closed on
#' exit, even if drawing fails.
#'
#' @param x An object of class \code{"mwperm"}.
#' @param file Output path; the extension selects the graphics device:
#'   \code{.pdf} (vector), \code{.png}, \code{.tiff}/\code{.tif}, or
#'   \code{.jpeg}/\code{.jpg}.
#' @param width \code{"single"} (3.5 in, a one-column journal figure),
#'   \code{"double"} (7 in, full text width), or a numeric width in inches.
#' @param height Height in inches; defaults to \code{0.65 * width}
#'   (\code{0.9 * width} for \code{type = "all"}).
#' @param type Figure type, passed to \code{\link{plot.mwperm}}.
#' @param res Raster resolution in dpi (default 300; ignored for pdf).
#' @param pointsize Base point size; defaults to 8 for single-column widths
#'   and 10 otherwise, so labels stay readable at print size.
#' @param ... Passed to \code{\link{plot.mwperm}} (style overrides, titles).
#' @return The file path, invisibly.
#' @seealso \code{\link{plot.mwperm}}.
#' @examples
#' data(trade_dyadic)
#' fit <- with(trade_dyadic,
#'             mwperm_dyadic(y = log_trade, d = log_dist,
#'                           x = cbind(log_gdp_i, log_gdp_j),
#'                           row = importer, col = exporter, seed = 1))
#' f <- file.path(tempdir(), "distance-elasticity.png")
#' mwperm_save(fit, f, width = "single")  # 3.5 in x 2.3 in at 300 dpi
#' mwperm_save(fit, file.path(tempdir(), "fig.pdf"),
#'             width = "double", type = "stability")
#' @export
mwperm_save <- function(x, file, width = c("single", "double"), height = NULL,
                        type = "auto", res = 300, pointsize = NULL, ...) {
  if (!inherits(x, "mwperm"))
    stop("`x` must be an object of class \"mwperm\".", call. = FALSE)
  if (!(is.character(file) && length(file) == 1L && nzchar(file)))
    stop("`file` must be a single file path.", call. = FALSE)
  if (is.character(width)) {
    width <- switch(match.arg(width), single = 3.5, double = 7)
  } else if (!(is.numeric(width) && length(width) == 1L && is.finite(width) &&
               width > 0)) {
    stop("`width` must be \"single\", \"double\", or a positive width in inches.",
         call. = FALSE)
  }
  if (is.null(height)) {
    height <- if (identical(type, "all")) 0.9 * width else 0.65 * width
  } else if (!(is.numeric(height) && length(height) == 1L &&
               is.finite(height) && height > 0)) {
    stop("`height` must be a positive height in inches.", call. = FALSE)
  }
  if (!(is.numeric(res) && length(res) == 1L && is.finite(res) && res >= 72))
    stop("`res` must be a single resolution >= 72 dpi.", call. = FALSE)
  if (is.null(pointsize)) pointsize <- if (width <= 4.5) 8 else 10
  ext <- tolower(sub(".*\\.", "", basename(file)))
  if (identical(ext, tolower(basename(file))))
    stop("`file` needs an extension: .pdf, .png, .tiff or .jpeg.", call. = FALSE)
  switch(ext,
         pdf  = grDevices::pdf(file, width = width, height = height,
                               pointsize = pointsize),
         png  = grDevices::png(file, width = width, height = height,
                               units = "in", res = res, pointsize = pointsize),
         tif  = ,
         tiff = grDevices::tiff(file, width = width, height = height,
                                units = "in", res = res, pointsize = pointsize,
                                compression = "lzw"),
         jpg  = ,
         jpeg = grDevices::jpeg(file, width = width, height = height,
                                units = "in", res = res, pointsize = pointsize,
                                quality = 95),
         stop(sprintf("Unsupported extension \".%s\"; use .pdf, .png, .tiff or .jpeg.",
                      ext), call. = FALSE))
  on.exit(grDevices::dev.off(), add = TRUE)
  plot(x, type = type, ...)
  invisible(file)
}
