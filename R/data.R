## ============================================================================
## R/data.R -- documentation of the bundled synthetic data sets
##
## Purpose. roxygen for trade_dyadic (a complete 40 x 40 dyadic array, two-way
##   random effects, GTW Eq. 8) and trade_panel (a 22 x 22 x 6 panel with a
##   common time trend, condition InvB). No code; data-raw/make_data.R makes
##   both (the JSS draft, Section 4).
## Pipeline. Not part of it: example inputs for every stage.
## ============================================================================
##
## Documentation for the bundled synthetic datasets. Both are generated
## reproducibly by data-raw/make_data.R (seeded; no external inputs).

#' Synthetic dyadic trade (gravity) data
#'
#' A simulated cross-section of bilateral trade among 40 countries, generated
#' from a gravity equation with a two-way (importer + exporter) clustered
#' error structure eps_ij = eta_i + xi_j + u_ij and heavy-tailed (scaled t_4)
#' idiosyncratic shocks. Designed to exercise [mwperm_dyadic()]; dropping the
#' self-trade diagonal yields an incomplete array suitable for
#' [mwperm_missing()].
#'
#' The data are synthetic and contain no real trade statistics. They are
#' generated reproducibly by the script in `data-raw/make_data.R`.
#'
#' @format A data frame with 1600 rows (a complete 40 x 40 array, including
#'   self-trade where `importer == exporter`) and 8 variables:
#' - `importer`: Factor with 40 levels (ISO3 codes); the row cluster.
#' - `exporter`: Factor with 40 levels (ISO3 codes); the column cluster.
#' - `log_trade`: Numeric outcome: log bilateral trade flow.
#' - `log_dist`: Numeric dyad-level covariate: log geographic distance between
#'   the pair (true coefficient -1.0).
#' - `log_gdp_i`: Numeric node-level covariate: log GDP of the importer (true
#'   coefficient 0.7); constant across exporters within an importer.
#' - `log_gdp_j`: Numeric node-level covariate: log GDP of the exporter (true
#'   coefficient 0.7).
#' - `border`: Integer 0/1 dyad-level covariate: shared land border (true
#'   coefficient 0.4).
#' - `placebo`: Numeric dyad-level covariate independent of the outcome (true
#'   coefficient 0); useful for checking the null rejection rate.
#'
#' The attribute `"true_coef"` stores the data-generating coefficients
#' (intercept 2.0, `log_dist` -1.0, `log_gdp_i` and `log_gdp_j` 0.7, `border`
#' 0.4, `placebo` 0).
#' @source Seeded simulation; the generator script is `data-raw/make_data.R`
#'   in the package's source repository (not installed with the package).
#' @seealso [mwperm_dyadic()], [mwperm_missing()], [mwperm()].
#' @examples
#' data(trade_dyadic)
#' str(trade_dyadic)
#' attr(trade_dyadic, "true_coef")
"trade_dyadic"

#' Synthetic panel (longitudinal) dyadic trade data
#'
#' A simulated panel of bilateral trade among 22 countries over six years. The
#' error structure adds an *arbitrary* (irregular) common time trend to the
#' two-way clustered components, eps_ijt = eta_i + xi_j + zeta_t + u_ijt, with
#' heavy-tailed (scaled t_4) idiosyncratic shocks. This is the setting of
#' condition InvB and [mwperm_panel()]: errors are exchangeable across
#' importer and exporter within each year but not over time.
#'
#' The data are synthetic and contain no real trade statistics. They are
#' generated reproducibly by the script in `data-raw/make_data.R`.
#'
#' @format A data frame with 2904 rows (a complete balanced 22 x 22 x 6 array)
#'   and 8 variables:
#' - `importer`: Factor with 22 levels (ISO3 codes); the row cluster.
#' - `exporter`: Factor with 22 levels (ISO3 codes); the column cluster.
#' - `year`: Integer time index, 2015--2020.
#' - `log_trade`: Numeric outcome: log bilateral trade flow.
#' - `fta`: Integer 0/1 dyad-by-time policy covariate: an in-force free trade
#'   agreement for the pair in that year (true coefficient 0.5). Roll-out
#'   timing varies randomly across pairs.
#' - `log_gdp_i`: Numeric time-varying node covariate: log GDP of the importer
#'   (true coefficient 0.6).
#' - `log_gdp_j`: Numeric time-varying node covariate: log GDP of the exporter
#'   (true coefficient 0.6).
#' - `placebo`: Numeric dyad-by-time covariate independent of the outcome
#'   (true coefficient 0).
#'
#' The attribute `"true_coef"` stores the data-generating coefficients
#' (intercept 1.5, `fta` 0.5, `log_gdp_i` and `log_gdp_j` 0.6, `placebo` 0).
#' @source Seeded simulation; the generator script is `data-raw/make_data.R`
#'   in the package's source repository (not installed with the package).
#' @seealso [mwperm_panel()], [mwperm()].
#' @examples
#' data(trade_panel)
#' str(trade_panel)
#' attr(trade_panel, "true_coef")
"trade_panel"
