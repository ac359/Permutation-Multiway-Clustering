## Documentation for the bundled synthetic datasets. Both are generated
## reproducibly by data-raw/make_data.R (seeded; no external inputs).

#' Synthetic dyadic trade (gravity) data
#'
#' A simulated cross-section of bilateral trade among 40 countries, generated
#' from a gravity equation with a two-way (importer + exporter) clustered error
#' structure \eqn{\varepsilon_{ij} = \eta_i + \xi_j + u_{ij}} and heavy-tailed
#' (scaled \eqn{t_4}) idiosyncratic shocks. Designed to exercise
#' \code{\link{mwperm_dyadic}}; dropping the self-trade diagonal yields an
#' incomplete array suitable for \code{\link{mwperm_missing}}.
#'
#' The data are synthetic and contain no real trade statistics. They are
#' generated reproducibly by the script in \code{data-raw/make_data.R}.
#'
#' @format A data frame with 1600 rows (a complete 40 x 40 array, including
#' self-trade where \code{importer == exporter}) and 8 variables:
#' \describe{
#'   \item{importer}{Factor with 40 levels (ISO3 codes); the row cluster.}
#'   \item{exporter}{Factor with 40 levels (ISO3 codes); the column cluster.}
#'   \item{log_trade}{Numeric outcome: log bilateral trade flow.}
#'   \item{log_dist}{Numeric dyad-level covariate: log geographic distance
#'     between the pair (true coefficient \eqn{-1.0}).}
#'   \item{log_gdp_i}{Numeric node-level covariate: log GDP of the importer
#'     (true coefficient \eqn{0.7}); constant across exporters within an
#'     importer.}
#'   \item{log_gdp_j}{Numeric node-level covariate: log GDP of the exporter
#'     (true coefficient \eqn{0.7}).}
#'   \item{border}{Integer 0/1 dyad-level covariate: shared land border (true
#'     coefficient \eqn{0.4}).}
#'   \item{placebo}{Numeric dyad-level covariate independent of the outcome
#'     (true coefficient \eqn{0}); useful for checking the null rejection rate.}
#' }
#' The attribute \code{"true_coef"} stores the data-generating coefficients
#' (intercept \eqn{2.0}, \code{log_dist} \eqn{-1.0}, \code{log_gdp_i} and
#' \code{log_gdp_j} \eqn{0.7}, \code{border} \eqn{0.4}, \code{placebo}
#' \eqn{0}).
#' @source Seeded simulation; the generator script is
#'   \code{data-raw/make_data.R} in the package's source repository (not
#'   installed with the package).
#' @examples
#' data(trade_dyadic)
#' str(trade_dyadic)
#' attr(trade_dyadic, "true_coef")
"trade_dyadic"

#' Synthetic panel (longitudinal) dyadic trade data
#'
#' A simulated panel of bilateral trade among 22 countries over six years. The
#' error structure adds an \emph{arbitrary} (irregular) common time trend to the
#' two-way clustered components, \eqn{\varepsilon_{ijt} = \eta_i + \xi_j +
#' \zeta_t + u_{ijt}}, with heavy-tailed (scaled \eqn{t_4}) idiosyncratic
#' shocks. This is the setting of condition InvB and
#' \code{\link{mwperm_panel}}: errors are exchangeable across importer and
#' exporter within each year but not over time.
#'
#' The data are synthetic and contain no real trade statistics. They are
#' generated reproducibly by the script in \code{data-raw/make_data.R}.
#'
#' @format A data frame with 2904 rows (a complete balanced 22 x 22 x 6 array)
#' and 8 variables:
#' \describe{
#'   \item{importer}{Factor with 22 levels (ISO3 codes); the row cluster.}
#'   \item{exporter}{Factor with 22 levels (ISO3 codes); the column cluster.}
#'   \item{year}{Integer time index, 2015--2020.}
#'   \item{log_trade}{Numeric outcome: log bilateral trade flow.}
#'   \item{fta}{Integer 0/1 dyad-by-time policy covariate: an in-force free
#'     trade agreement for the pair in that year (true coefficient \eqn{0.5}).
#'     Roll-out timing varies randomly across pairs.}
#'   \item{log_gdp_i}{Numeric time-varying node covariate: log GDP of the
#'     importer (true coefficient \eqn{0.6}).}
#'   \item{log_gdp_j}{Numeric time-varying node covariate: log GDP of the
#'     exporter (true coefficient \eqn{0.6}).}
#'   \item{placebo}{Numeric dyad-by-time covariate independent of the outcome
#'     (true coefficient \eqn{0}).}
#' }
#' The attribute \code{"true_coef"} stores the data-generating coefficients
#' (intercept \eqn{1.5}, \code{fta} \eqn{0.5}, \code{log_gdp_i} and
#' \code{log_gdp_j} \eqn{0.6}, \code{placebo} \eqn{0}).
#' @source Seeded simulation; the generator script is
#'   \code{data-raw/make_data.R} in the package's source repository (not
#'   installed with the package).
#' @examples
#' data(trade_panel)
#' str(trade_panel)
#' attr(trade_panel, "true_coef")
"trade_panel"
