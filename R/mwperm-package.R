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
#' balance, attainable resolution) without running anything. The
#' design-specific tests -- which [mwperm()] calls and which remain fully
#' supported for direct use -- are [mwperm_dyadic()] (two-way / dyadic
#' clustering), [mwperm_threeway()] (three-way clustering), [mwperm_panel()]
#' (panels with an arbitrary time effect), [mwperm_layout()] (replicated
#' two-way layouts), [mwperm_irregular()] (irregular layouts, Section 6.4:
#' exchangeable replicates with a covariate constant within cells),
#' [mwperm_missing()] (incomplete arrays, via fully observed bicliques) and
#' [mwperm_dyadic_het()] (dyadic clustering under heteroskedasticity). All
#' return an object of class `"mwperm"` with [print.mwperm()],
#' [summary.mwperm()], [confint.mwperm()] and [plot.mwperm()] methods.
#'
#' Every test is the same procedure -- partial out the nuisance design
#' against its transformed copy, minorize, count -- run under a group of
#' transformations the errors are assumed invariant to, and the package has
#' two group constructions with two different assumptions:
#' - **Permutations** ([build_perm_set()], the block-cyclic group of
#'   Algorithm 1): the error array must be *exchangeable* under relabelling
#'   of the clusters, conditional on the covariates (Assumption 1). Used by
#'   every test except the last. Tolerates any error distribution and any
#'   dependence that is symmetric in the cluster labels (additive cluster
#'   effects included); does not tolerate an error variance, or any other
#'   feature of the error law, that depends on the covariates.
#' - **Sign flips** ([build_flip_set()], the group of joint row-and-column
#'   sign changes, order `2^(n_flip - 1)`): the error array must be
#'   *jointly symmetric* under those sign changes. Used by
#'   [mwperm_dyadic_het()]. Tolerates arbitrary heteroskedasticity for
#'   errors independent across cells; does not tolerate skewness, and does
#'   **not** tolerate additive cluster effects `eta_i + xi_j`, which change
#'   the joint law under a sign flip (the test then over-rejects).
#'
#' Neither assumption implies the other, and the sign-flip test is less
#' powerful when exchangeability does hold, so the choice is a judgement
#' about the errors that the data cannot make for you: [mwperm()] never
#' selects the sign-flip test automatically.
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
#' @seealso [mwperm()], [mwperm_check()], [mwperm_dyadic()],
#'   [mwperm_dyadic_het()], [mwperm_panel()], [build_perm_set()],
#'   [build_flip_set()].
#'
#' @importFrom stats median lm.fit model.matrix sd setNames coef nobs
#' @importFrom graphics arrows axis hist legend mtext par plot.new plot.window
#'   points polygon rect segments strwidth text title
#' @importFrom grDevices adjustcolor chull dev.off jpeg pdf png tiff
#' @keywords internal
"_PACKAGE"
