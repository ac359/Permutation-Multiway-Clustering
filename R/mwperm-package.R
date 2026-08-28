#' mwperm: Invariant Permutation Tests for Multi-Way Clustered and Panel
#' Regression
#'
#' Finite-sample valid tests and confidence intervals for regression
#' coefficients under multi-way (e.g. dyadic) clustering, including panel and
#' missing-data designs, implementing the invariant permutation test of Guo,
#' Toulis and Wang (2026).
#'
#' The recommended entry points are [mwperm()], which detects the clustering
#' design of the data and dispatches to the matching test, and
#' [mwperm_check()], which prints the diagnosis (detected design, roles,
#' balance, attainable resolution) without running anything. The key object is
#' the block-cyclic permutation group built by [build_perm_set()]. The
#' design-specific tests -- which [mwperm()] calls and which remain fully
#' supported for direct use -- are [mwperm_dyadic()] (two-way / dyadic
#' clustering), [mwperm_threeway()] (three-way clustering), [mwperm_panel()]
#' (panels with an arbitrary time effect), [mwperm_layout()] (replicated
#' two-way layouts), [mwperm_irregular()] (irregular layouts, Section 6.4:
#' repeats that are time periods, or a covariate constant within cells) and
#' [mwperm_missing()] (incomplete arrays, via fully observed bicliques). All
#' return an object of class `"mwperm"` with [print.mwperm()],
#' [summary.mwperm()], [confint.mwperm()] and [plot.mwperm()] methods.
#'
#' Two synthetic data sets, [trade_dyadic] and [trade_panel], illustrate the
#' dyadic and panel work flows.
#'
#' @section Reproducibility: A `seed` determines the result completely: it
#'   draws the permutation group, and every step after that is deterministic.
#'   For a fixed seed, platform and BLAS, repeated runs return the same
#'   p-values, estimates and interval endpoints bit for bit. Under a
#'   *different* BLAS the underlying matrix products are blocked and
#'   reassociated differently, so endpoints may differ in their final
#'   decimals. The p-value is unaffected in any practical sense: it lives on
#'   the grid `{1, ..., K+1}/(K+1)`, whose resolution 1/(K+1) is coarser than
#'   such perturbations by many orders of magnitude.
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data. arXiv:2601.08610.
#'
#' @seealso [mwperm()], [mwperm_check()], [mwperm_dyadic()], [mwperm_panel()],
#'   [build_perm_set()].
#'
#' @importFrom stats median lm.fit model.matrix sd setNames coef nobs
#' @importFrom graphics arrows axis hist legend mtext par plot.new plot.window
#'   points polygon rect segments strwidth text title
#' @importFrom grDevices adjustcolor chull dev.off jpeg pdf png tiff
#' @keywords internal
"_PACKAGE"
