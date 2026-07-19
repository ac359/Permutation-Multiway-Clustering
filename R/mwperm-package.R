#' mwperm: Invariant Permutation Tests for Multi-Way Clustered and Panel
#'   Regression
#'
#' Finite-sample valid tests and confidence intervals for regression
#' coefficients under multi-way (e.g. dyadic) clustering, including panel and
#' missing-data designs, implementing the invariant permutation test of Guo,
#' Toulis and Wang (2026).
#'
#' The recommended entry points are \code{\link{mwperm}}, which detects the
#' clustering design of the data and dispatches to the matching test, and
#' \code{\link{mwperm_check}}, which prints the diagnosis (detected design,
#' roles, balance, attainable resolution) without running anything. The key
#' object is the block-cyclic permutation group built by
#' \code{\link{build_perm_set}}. The design-specific tests -- which
#' \code{\link{mwperm}} calls and which remain fully supported for direct
#' use -- are \code{\link{mwperm_dyadic}} (two-way / dyadic clustering),
#' \code{\link{mwperm_threeway}} (three-way clustering),
#' \code{\link{mwperm_panel}} (panels with an arbitrary time effect),
#' \code{\link{mwperm_layout}} (replicated two-way layouts) and
#' \code{\link{mwperm_missing}} (incomplete arrays, via fully observed
#' bicliques). All return an object of class \code{"mwperm"} with
#' \code{\link{print.mwperm}}, \code{\link{summary.mwperm}},
#' \code{\link{confint.mwperm}} and \code{\link{plot.mwperm}} methods.
#'
#' Two synthetic data sets, \code{\link{trade_dyadic}} and
#' \code{\link{trade_panel}}, illustrate the dyadic and panel work flows.
#'
#' @references
#' Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference under
#' multi-way clustering and missing data. arXiv:2601.08610.
#'
#' @seealso \code{\link{mwperm}}, \code{\link{mwperm_check}},
#'   \code{\link{mwperm_dyadic}}, \code{\link{mwperm_panel}},
#'   \code{\link{build_perm_set}}.
#'
#' @importFrom stats median lm.fit model.matrix sd setNames coef nobs
#' @importFrom graphics arrows axis hist legend mtext par plot.new
#'   plot.window points polygon rect segments strwidth text title
#' @importFrom grDevices adjustcolor chull dev.off jpeg pdf png tiff
#' @keywords internal
"_PACKAGE"
