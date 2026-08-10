## Formula interface: a thin assembly layer over mwperm(). The statistical
## core is untouched -- the wrapper only builds y/d/x from formula algebra
## and forwards; identity with the data interface is pinned by
## tests/test-formula.R.

#' Formula interface to the invariant permutation test
#'
#' Fits \code{\link{mwperm}} from a formula \code{y ~ d | x} (the nuisance
#' part is optional: \code{y ~ d}): the part between \code{~} and \code{|}
#' gives the covariate(s) of interest, the part after \code{|} the nuisance
#' covariates. Both parts are standard formula algebra evaluated by
#' \code{\link[stats]{model.matrix}} in \code{data}, so transformed terms
#' (\code{log(z)}, \code{I(z^2)}, factors) work; intercept columns are
#' stripped (the engine always adds its own). The result is
#' \emph{identical} to calling \code{\link{mwperm}} -- or the dispatched
#' design-specific function -- with the same inputs and seed; this wrapper
#' only assembles the arguments.
#'
#' @param formula A two-sided formula, \code{y ~ d} or \code{y ~ d | x}.
#' @param data A data frame in which the formula (and character
#'   \code{index}/\code{time}/\code{rep}) are evaluated.
#' @param index The clustering dimensions (2 or 3): a character vector of
#'   column names in \code{data}, or a data frame / named list of vectors.
#' @param time,rep Optional role tags (a column name in \code{data}, or a
#'   vector); see \code{\link{mwperm_check}}.
#' @param ... Passed on to \code{\link{mwperm}} (\code{design}, \code{K},
#'   \code{alpha}, \code{beta_null}, \code{conf_int}, \code{n_reps},
#'   \code{seed}, \code{verbose}, ...).
#' @return The \code{"mwperm"} object of the dispatched test; see
#'   \code{\link{mwperm}} for the fields and their provenance.
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation
#'   inference under multi-way clustering and missing data. arXiv:2601.08610.
#' @seealso \code{\link{mwperm}}; \code{\link{coef.mwperm}} /
#'   \code{\link{nobs.mwperm}} for accessors.
#' @examples
#' data(trade_dyadic)
#' fit <- mwperm_formula(log_trade ~ log_dist | log_gdp_i + log_gdp_j,
#'                       data = trade_dyadic,
#'                       index = c("importer", "exporter"),
#'                       n_reps = 3, seed = 1, verbose = FALSE)
#' fit
#' coef(fit)
#' nobs(fit)
#' @export
mwperm_formula <- function(formula, data, index, time = NULL, rep = NULL, ...) {
  if (!(inherits(formula, "formula") && length(formula) == 3L))
    stop("`formula` must be two-sided: y ~ d or y ~ d | x.", call. = FALSE)
  rhs <- formula[[3L]]
  if (is.call(rhs) && identical(rhs[[1L]], as.name("|"))) {
    d_part <- rhs[[2L]]
    x_part <- rhs[[3L]]
  } else {
    d_part <- rhs
    x_part <- NULL
  }
  mm <- function(part) {
    f <- stats::as.formula(call("~", part), env = environment(formula))
    m <- stats::model.matrix(f, data = data)
    m[, colnames(m) != "(Intercept)", drop = FALSE]
  }
  y <- eval(formula[[2L]], data, environment(formula))
  d <- mm(d_part)
  x <- if (is.null(x_part)) NULL else mm(x_part)
  idx <- if (is.character(index)) data[index] else as.data.frame(index)
  tv  <- if (is.character(time) && length(time) == 1L) data[[time]] else time
  rv  <- if (is.character(rep)  && length(rep)  == 1L) data[[rep]]  else rep
  mwperm(y = y, d = d, x = x, index = idx, time = tv, rep = rv, ...)
}

#' Accessors for mwperm fits
#'
#' \code{coef()} returns the \emph{OLS} point estimate(s), named by
#' coefficient -- the same provenance as the printed \code{"OLS estimate"};
#' the inferential quantity is the IPT confidence set, see
#' \code{\link{confint.mwperm}}. \code{nobs()} returns the number of
#' observations the fit used (for \code{\link{mwperm_missing}}, the cells
#' inside the selected blocks).
#'
#' @param object An object of class \code{"mwperm"}.
#' @param ... Ignored.
#' @return \code{coef()}: a named numeric vector; \code{nobs()}: an integer.
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation
#'   inference under multi-way clustering and missing data. arXiv:2601.08610.
#' @seealso \code{\link{confint.mwperm}}, \code{\link{summary.mwperm}}.
#' @examples
#' data(trade_dyadic)
#' fit <- mwperm_formula(log_trade ~ log_dist | log_gdp_i + log_gdp_j,
#'                       data = trade_dyadic,
#'                       index = c("importer", "exporter"),
#'                       n_reps = 2, seed = 1, verbose = FALSE)
#' coef(fit)
#' nobs(fit)
#' @name mwperm-accessors
NULL

#' @rdname mwperm-accessors
#' @export
coef.mwperm <- function(object, ...)
  stats::setNames(as.numeric(object$estimate), object$d_names)

#' @rdname mwperm-accessors
#' @export
nobs.mwperm <- function(object, ...) object$n_obs
