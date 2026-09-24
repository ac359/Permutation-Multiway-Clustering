## Unified entry point: automatic design detection (mwperm_check) and
## dispatch (mwperm). A thin, additive layer over the eight mwperm_* front
## ends -- it changes nothing about how any test is computed.
##
## Detection policy. Forks that are visible in the data's *structure*
## (index count, repeated cells, completeness) are resolved silently; forks
## that hinge on an EXCHANGEABILITY ASSUMPTION the data cannot reveal are
## announced, default to the choice that stays valid under the widest set of
## data-generating processes, and carry override instructions. The asymmetry
## is deliberate: structure is checkable from the data, exchangeability is
## not, so guessing wrong about structure fails loudly while guessing wrong
## about exchangeability fails silently and invalidates the test.
##   * 3 complete crossed indices: panel vs threeway is such a fork. Running
##     threeway on a panel (time-autocorrelated errors) is INVALID (size
##     distortion; in this package's own Monte Carlo, 0.88 rejection at
##     alpha = .2 on a trending panel null against the panel test's 0.197);
##     running panel on genuinely three-way
##     exchangeable data is merely less powerful (InvA implies InvB). Hence
##     panel is the default; threeway only on explicit design = "threeway".
##   * 2 indices with repeated cells: layout vs a panel whose time index was
##     not passed. Layout assumes the within-cell replicates are exchangeable;
##     if they are a time series that is false. Default layout + warning.
##   * Heteroskedasticity is a third such assumption, and it has NO
##     structural signature at all, so the sign-flip test (mwperm_dyadic_het,
##     design = "dyadic_het") is never detected -- only offered, in one line
##     of the printed diagnosis on a complete dyadic array. That line lives in
##     `alternatives`, not `notes`, so it never reaches a fitted object.

## ---- small helpers ----------------------------------------------------------

#' Resolve the `index` argument to a named list of equal-length vectors.
#' @keywords internal
#' @noRd
.resolve_index <- function(index, data) {
  if (is.character(index) && is.null(dim(index))) {
    ## character vector of column names against `data`
    if (is.null(data))
      stop(paste0("`index` is a character vector of column names but ",
                  "`data` was not supplied."),
           call. = FALSE)
    bad <- setdiff(index, names(data))
    if (length(bad))
      stop(sprintf("Index column(s) %s not found in `data`.",
                   paste0("'", bad, "'", collapse = ", ")), call. = FALSE)
    out <- as.list(data[index])
    names(out) <- index
    return(out)
  }
  if (is.data.frame(index) || is.list(index)) {
    out <- as.list(index)
    nm <- names(out)
    if (is.null(nm)) nm <- rep("", length(out))
    nm[!nzchar(nm)] <- paste0("index", which(!nzchar(nm)))
    names(out) <- nm
    return(out)
  }
  stop("`index` must be a data frame, a (named) list of vectors, or a ",
       "character vector of column names in `data`.", call. = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Is a dimension name time-like? (case-insensitive vocabulary)
#' @keywords internal
#' @noRd
.timelike_name <- function(nm) {
  if (is.null(nm) || !nzchar(nm)) return(FALSE)
  vocab <- c("year", "yr", "time", "t", "date", "period", "wave", "month",
             "quarter", "day", "decade", "week", "season", "annum")
  tolower(nm) %in% c(vocab, paste0(vocab, "s"))
}

#' Are a dimension's values time-like? (temporal class, or regularly spaced
#' numeric with no more levels than the cross-sectional dimensions)
#' @keywords internal
#' @noRd
.timelike_values <- function(v, n_levels, max_levels) {
  if (inherits(v, c("Date", "POSIXct", "POSIXlt"))) return(TRUE)
  if (!is.numeric(v)) return(FALSE)
  u <- sort(unique(as.numeric(v)))
  if (length(u) < 2L || length(u) > max_levels) return(FALSE)
  d <- diff(u)
  isTRUE(all(abs(d - d[1L]) < 1e-8))             # regularly spaced
}

#' Strongly time-like values: temporal class, or regularly spaced numeric with
#' strictly fewer levels than every other dimension (the typical few-periods
#' panel shape). Consecutive integer cluster ids satisfy the weak rule above,
#' so this stricter form is what decides whether a NAME-based time assignment
#' is corroborated by the values. Used only to gate warnings -- never to
#' assign the time role itself (assignment behaviour is frozen).
#' @keywords internal
#' @noRd
.timelike_strong <- function(v, n_levels, min_other) {
  if (inherits(v, c("Date", "POSIXct", "POSIXlt"))) return(TRUE)
  if (!is.numeric(v) || n_levels >= min_other) return(FALSE)
  u <- sort(unique(as.numeric(v)))
  length(u) >= 2L && isTRUE(all(abs(diff(u) - diff(u)[1L]) < 1e-8))
}

## ---- mwperm_check -----------------------------------------------------------

#' Diagnose a multi-way clustered dataset and choose the appropriate design
#'
#' Inspects the clustering structure of a dataset -- number of index
#' dimensions, repeated cells, completeness/balance -- and reports which
#' `mwperm_*` test applies, without running any permutations or fitting
#' anything. [mwperm()] uses it for automatic dispatch; call it directly to
#' see the diagnosis.
#'
#' Structural forks (dyadic vs missing, panel vs incomplete panel, and
#' layout-by-replication) are resolved silently from the data. Two forks
#' depend on an *exchangeability assumption the data cannot reveal* and are
#' therefore announced with override instructions, defaulting to the choice
#' that remains valid under the widest set of error processes:
#' - **panel vs three-way** (complete balanced 3-index arrays): running
#'   [mwperm_threeway()] on a panel whose errors are dependent over time is
#'   *invalid* (size distortion), while running [mwperm_panel()] on genuinely
#'   three-way exchangeable data is valid, only less powerful. The default is
#'   therefore `panel`. The time role is assigned by, in order: an explicit
#'   `time =` tag; a time-like *name* (case-insensitive vocabulary: year, yr,
#'   time, t, date, period, wave, month, quarter, day, decade, week, season,
#'   annum, and their plurals); time-like *values* (temporal class, or
#'   regularly spaced numeric). A name-based assignment that the values do not
#'   corroborate (temporal class, or regularly spaced with strictly fewer
#'   levels than every other dimension) carries a **warning** -- a column
#'   merely *named* like time may be a cluster, and permuting the true time
#'   dimension over-rejects badly. An ambiguous case defaults to holding the
#'   third index fixed, with a warning. Forcing `design = "threeway"` when an
#'   index looks time-like also warns; force it only when all three dimensions
#'   are genuinely exchangeable.
#' - **layout vs suppressed panel** (2 indices with repeated cells): repeats
#'   are treated as within-cell replication ([mwperm_layout()]), which assumes
#'   the replicates are exchangeable within cells -- if they are really a time
#'   series, pass the time variable via `time =` to get the panel test
#'   instead. A notice is attached.
#'
#' A third assumption the data cannot reveal is *heteroskedasticity*: every
#' permutation design needs the errors exchangeable *given* the covariates,
#' which an error variance that depends on the covariates violates. Nothing
#' in the clustering structure shows this, so the diagnosis never selects
#' the sign-flip test ([mwperm_dyadic_het()], valid under arbitrary
#' heteroskedasticity for errors that are independent across cells and
#' symmetric, but not under additive cluster effects); on a dyadic array,
#' complete or not, it prints one line offering `design = "dyadic_het"`, and
#' the choice is yours.
#'
#' @param index The clustering dimensions (2 or 3): a data frame, a named list
#'   of vectors, or a character vector of column names resolved against
#'   `data`.
#' @param y,d Optional outcome and covariate(s) of interest; only used for
#'   extra diagnostics (e.g. the layout no-power warning when `d` is constant
#'   within every cell), never for fitting.
#' @param data Optional data frame against which character `index`, `time` and
#'   `rep` entries are resolved.
#' @param time Optional explicit time dimension: a vector, or the name of a
#'   column of `data` (or of one of the `index` columns). Forces the panel
#'   interpretation of that dimension: a complete `(row, col, time)` array
#'   runs [mwperm_panel()], an incomplete one [mwperm_panel_missing()] (to
#'   which [mwperm()] also forwards `L0 =`, the number of periods to keep
#'   when some pairs miss periods). With `design = "irregular"` it runs
#'   [mwperm_panel_missing()] too: periods are never given the random
#'   per-cell trim, which is exact only for exchangeable replicates.
#' @param rep Optional explicit replication identifier (vector or column
#'   name): declares within-cell replication and forces the layout design.
#' @param design Force a design instead of auto-detecting (the structure is
#'   still validated against it). `"dyadic_het"` -- the sign-flip test of
#'   [mwperm_dyadic_het()] -- is *opt-in only*: it is never detected, because
#'   heteroskedasticity leaves no trace in the clustering structure, and it
#'   is validated exactly as `"dyadic"` (two indices, one observation per
#'   cell, complete array).
#' @param alpha,aggregate The test level and cross-repetition rule the fit
#'   will use (the defaults of every front end). They decide the resolution
#'   verdict: the smallest reportable p-value is `1/(K+1)` under `"median"`
#'   and `2/(K+1)` under `"median2"`, and a `(1 - alpha)` confidence set is
#'   attainable only when that floor is at most `alpha`. [mwperm()] passes
#'   its own `alpha` and `aggregate` through, so the diagnosis it prints
#'   describes the fit it runs.
#'
#' @return An object of class `"mwperm_design"`: a list with fields `design`
#'   (the chosen design), `roles` (which index plays row/col/id1..3/time/rep),
#'   `dims` (levels per dimension), `n_obs`, `cells` (observed/expected),
#'   `balance`, `K_default`, `alpha`, `aggregate`, `p_floor` (the smallest
#'   reportable p-value at the default K), `levels_needed` (the smallest
#'   permuted dimension a `(1 - alpha)` set requires) and `resolution_ok`
#'   (whether that set is attainable), `call_str` (the downstream call),
#'   `reason` (one-line explanation), `warnings`/`notes` (the
#'   assumption-fork notices etc.) and `alternatives` (one-line pointers to
#'   designs the data cannot select for you -- on a complete dyadic array,
#'   `design = "dyadic_het"`). For `design = "dyadic_het"` the group order
#'   is `2^(n_flip - 1)`, so `K_default` is `NA`, the extra field
#'   `n_flip_default` carries the fit's default `n_flip` (the smallest whose
#'   p-value floor is at most `alpha` under the given `aggregate`, capped by
#'   the smaller dimension: 6 at `alpha = 0.05`), and `levels_needed` is the
#'   `n_flip` a `(1 - alpha)` set requires. Its `print` method lays this out
#'   as a short human diagnosis.
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data. arXiv:2601.08610.
#' @seealso [mwperm()] for one-call dispatch.
#' @examples
#' data(trade_dyadic)
#' mwperm_check(index = c("importer", "exporter"), data = trade_dyadic)
#' data(trade_panel)
#' mwperm_check(index = c("importer", "exporter", "year"), data = trade_panel)
#' @export
mwperm_check <- function(index, y = NULL, d = NULL, data = NULL,
                         time = NULL, rep = NULL,
                         design = c("auto", "dyadic", "threeway", "panel",
                                    "panel_missing", "layout", "missing",
                                    "irregular", "dyadic_het"),
                         alpha = 0.05, aggregate = c("median", "median2")) {
  design <- match.arg(design)
  aggregate <- match.arg(aggregate)
  if (!(is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
        alpha > 0 && alpha < 1))
    stop("`alpha` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  idx <- .resolve_index(index, data)
  N <- length(idx[[1L]])
  if (any(vapply(idx, length, integer(1)) != N))
    stop("All index dimensions must have the same length.", call. = FALSE)

  notes <- character(0)
  warns <- character(0)

  ## Explicit role tags. A character tag naming an index column CLAIMS that
  ## column (removes it from `idx` via <<-, so it no longer counts as a
  ## clustering dimension); a character tag naming a `data` column pulls that
  ## column in; a vector tag is used as-is. Returns a named one-element list
  ## (the name labels the role in roles/dims) or NULL.
  claim <- function(tag, what) {
    if (is.null(tag)) return(NULL)
    if (is.character(tag) && length(tag) == 1L) {
      if (tag %in% names(idx)) {
        v <- idx[[tag]]
        idx[[tag]] <<- NULL
        return(stats::setNames(list(v), tag))
      }
      if (!is.null(data) && tag %in% names(data))
        return(stats::setNames(list(data[[tag]]), tag))
      stop(sprintf(paste0("`%s = \"%s\"` matches neither an index column ",
                          "nor a column of `data`."),
                   what, tag), call. = FALSE)
    }
    if (length(tag) != N)
      stop(sprintf("`%s` must have the same length as the index columns (%d).",
                   what, N), call. = FALSE)
    stats::setNames(list(tag), what)
  }
  time_v <- claim(time, "time")
  rep_v <- claim(rep, "rep")

  ## drop single-level index columns (they carry no clustering)
  lvls <- vapply(idx, function(v) length(unique(v)), integer(1))
  if (any(lvls < 2L)) {
    dropped <- names(idx)[lvls < 2L]
    notes <- c(notes, sprintf(
      paste0("Index dimension(s) %s have a single level and were dropped ",
             "(no clustering)."),
      paste0("'", dropped, "'", collapse = ", ")))
    idx <- idx[lvls >= 2L]
  }
  C <- length(idx)
  if (C < 2L)
    stop("Fewer than 2 effective clustering dimensions remain. mwperm needs 2 ",
         "or 3 multi-level index dimensions (plus an optional time/rep role).",
         call. = FALSE)
  if (C > 3L || (C == 3L && (!is.null(time_v) || !is.null(rep_v))))
    stop("Too many clustering dimensions: pass 2 or 3 index columns ",
         "(a tagged `time` or `rep` counts as the third).", call. = FALSE)

  dims <- vapply(idx, function(v) length(unique(v)), integer(1))
  dense <- Map(.dense_id, idx, names(idx))
  cells2 <- .cell_code(cbind(dense[[1L]], dense[[2L]]))
  dup2 <- anyDuplicated(cells2) > 0L

  ## defaults filled per design below
  roles <- NULL
  chosen <- NULL
  reason <- NULL
  K_default <- NA_integer_
  balance <- NA_character_
  cells_obs <- NA_integer_
  cells_exp <- NA_integer_
  ## Sign-flip design only: the default n_flip the fit would use (its group
  ## order is 2^(n_flip - 1), so K_default is not the right number to show).
  n_flip_default <- NULL
  ## One-line pointers to a design the data cannot select for the user.
  alternatives <- character(0)

  finish_layout <- function(why, warn_txt = NULL) {
    cell <- .dense_id(interaction(dense[[1L]], dense[[2L]], drop = TRUE))
    sizes <- tabulate(cell)
    chosen <<- "layout"
    roles <<- list(row = names(idx)[1L], col = names(idx)[2L],
                   rep = if (!is.null(rep_v)) names(rep_v) else
                     "(within-cell order)")
    reason <<- why
    K_default <<- min(sizes) - 1L
    ## "4-4 replicates" read as a range when there is none; and the smallest
    ## cell is what sets K, which is the number the reader actually needs.
    balance <<- if (min(sizes) == max(sizes))
      sprintf("replicated (%d cells, %d replicates each)",
              length(sizes), min(sizes))
    else
      sprintf(paste0("replicated (%d cells, %d-%d replicates; the smallest ",
                     "cell sets K)"), length(sizes), min(sizes), max(sizes))
    cells_obs <<- length(sizes)
    cells_exp <<- prod(dims[1:2])
    if (!is.null(warn_txt)) warns <<- c(warns, warn_txt)
    ## no-power diagnostic when d is available
    if (!is.null(d)) {
      D <- as.matrix(d)
      wv <- tapply(seq_len(N), cell, function(ii)
        any(apply(D[ii, , drop = FALSE], 2L, function(v) diff(range(v)) > 0)))
      if (!any(unlist(wv)))
        warns <<- c(warns, paste0(
          "`d` is constant within every (row, col) cell: the within-cell ",
          "layout test (Section 6.3) would have NO power, because permuting ",
          "inside a cell cannot move a covariate that is constant there. Use ",
          "design = \"irregular\" with an `L0 =` threshold instead: that is ",
          "the Section 6.4 procedure, which permutes the CELLS across rows ",
          "and columns and is designed for exactly this case when the ",
          "repeats are exchangeable replicates; if they are periods, pass ",
          "them as `time =` instead -- see ?mwperm_irregular)."))
    }
  }

  finish_irregular <- function(why) {
    cell <- .dense_id(interaction(dense[[1L]], dense[[2L]], drop = TRUE))
    sizes <- tabulate(cell)
    chosen <<- "irregular"
    roles <<- list(row = names(idx)[1L], col = names(idx)[2L],
                   rep = if (!is.null(rep_v)) names(rep_v) else
                     "(within-cell order)")
    reason <<- why
    ## K is set by the biclique blocks found under the Section 6.4 mask
    ## M_ij = 1{ell_ij >= L0}, which depends on L0 -- not knowable here.
    K_default <<- NA_integer_
    balance <<- sprintf(paste0("irregular (%d cells, %d-%d observations ",
                               "each; L0 sets which cells are usable)"),
                        length(sizes), min(sizes), max(sizes))
    cells_obs <<- length(sizes)
    cells_exp <<- prod(dims[1:2])
    notes <<- c(notes, paste0(
      "The permutation-group order for the Section 6.4 design is set by the ",
      "biclique blocks found under the mask M_ij = 1{ell_ij >= L0}, so it ",
      "depends on `L0`; see find_bicliques() and ?mwperm_irregular. Its ",
      "random per-cell trim is exact only for exchangeable replicates within ",
      "a cell; if the within-cell index is a period, or anything else with ",
      "an effect shared across cells, it can over-reject even with a ",
      "cell-constant `d` -- use `time =` with `L0 =` instead ",
      "(mwperm_panel_missing(), the same L0 periods in every cell)."))
  }

  ## Structural gate shared by every panel path (auto, tagged, forced). A
  ## REPEATED (row, col, time) cell is always fatal: two observations would
  ## claim one slot of the "same (pi, sigma) in every period" map, and no
  ## panel design is defined. An INCOMPLETE array is not fatal -- since 0.4.0
  ## it is the incomplete-panel design, mwperm_panel_missing() -- so this only
  ## reports completeness and leaves the caller to route (auto) or refuse
  ## (a forced design = "panel").
  panel_cells <- function() {
    nT <- length(unique(time_v[[1L]]))
    tri <- .cell_code(cbind(dense[[1L]], dense[[2L]],
                            .dense_id(time_v[[1L]], names(time_v))))
    exp3 <- prod(dims[1:2]) * nT
    if (anyDuplicated(tri) > 0L)
      stop(sprintf(paste0(
        "A panel needs one observation per (row, col, time) cell, but some ",
        "(%s, %s, %s) cell repeats: %d observations for %s=%d x %s=%d x ",
        "%s=%d = %d cells. If the repeats are exchangeable replication ",
        "rather than time, pass them as `rep =` for a layout design; ",
        "otherwise aggregate to one observation per cell."),
        names(idx)[1L], names(idx)[2L], names(time_v), N,
        names(idx)[1L], dims[1L], names(idx)[2L], dims[2L],
        names(time_v), nT, exp3), call. = FALSE)
    list(expected = exp3, complete = N == exp3,
         dims_str = sprintf("%s=%d x %s=%d x %s=%d", names(idx)[1L], dims[1L],
                            names(idx)[2L], dims[2L], names(time_v), nT))
  }

  ## With no `time =` tag, a forced panel design takes the third index as time
  ## (the caller has already required C == 3).
  third_index_as_time <- function() {
    if (!is.null(time_v)) return(invisible())
    time_v <<- idx[3L]
    idx <<- idx[1:2]
    dims <<- dims[1:2]
    dense <<- dense[1:2]
  }

  ## The incomplete-panel design, mwperm_panel_missing(): reached by
  ## auto-detection (three untagged indices, or a tagged `time =`, on an
  ## incomplete array) and by forcing it. `announce` attaches the routing
  ## note on the auto routes only -- a forced choice needs no explanation,
  ## and the fit's own note already reports what the mask kept.
  finish_panel_missing <- function(why, complete, expected, announce) {
    chosen <<- "panel_missing"
    roles <<- list(row = names(idx)[1L], col = names(idx)[2L],
                   time = names(time_v))
    reason <<- why
    ## The group order comes from the biclique blocks, which are not searched
    ## until fit time, so K is unknown here (as for "missing").
    K_default <<- NA_integer_
    balance <<- if (complete) "complete" else "incomplete"
    cells_obs <<- N
    cells_exp <<- expected
    if (announce)
      notes <<- c(notes, paste0(
        "The array is incomplete, so the test restricts to (row, col) pairs ",
        "observed in EVERY period and to the fully observed blocks the ",
        "biclique search extracts from them; the group order follows those ",
        "blocks. If some pairs miss periods, `L0 =` keeps instead the L0 ",
        "periods jointly observed by the most pairs, the same ones in every ",
        "cell. See ?mwperm_panel_missing and find_bicliques()."))
  }

  if (design != "auto") {
    ## forced design: detection is skipped but structural feasibility is
    ## still validated (a forced design that cannot run should fail here,
    ## actionably, not deep inside the engine)
    chosen <- design
    reason <- "forced via design ="
    if (design %in% c("dyadic", "missing", "dyadic_het")) {
      if (C != 2L)
        stop(sprintf("design = \"%s\" needs exactly 2 index dimensions.",
                     design), call. = FALSE)
      if (dup2)
        stop(sprintf(paste0("design = \"%s\" requires one observation per ",
                            "cell, but some (row, col) cells repeat. Use ",
                            "layout (rep) or panel (time)."),
                             design), call. = FALSE)
      roles <- list(row = names(idx)[1L], col = names(idx)[2L])
      K_default <- min(dims[1:2]) - 1L
      cells_obs <- N
      cells_exp <- prod(dims[1:2])
      balance <- if (N == cells_exp) "complete"
                 else sprintf("incomplete (%d of %d cells)", N, cells_exp)
      ## A sign flip moves no observation, so "dyadic_het" runs on an
      ## incomplete array as it is (every observed cell, nothing discarded);
      ## only the permutation design needs the complete array.
      if (design == "dyadic" && N != cells_exp)
        stop(sprintf(paste0("design = \"%s\" requires a complete ",
                            "array but only ",
                            "%d of %d cells are observed. Use ",
                            "design = \"missing\" ",
                            "(Procedure 2, fully observed bicliques) instead."),
                     design, N, cells_exp), call. = FALSE)
      if (design == "missing") {
        K_default <- NA_integer_          # depends on the biclique blocks
      }
      if (design == "dyadic_het") {
        ## The sign-flip group's order is 2^(n_flip - 1) at the fit's default
        ## n_flip -- the smallest whose floor clears alpha under the fit's
        ## aggregation, capped by the smaller dimension: .default_n_flip(),
        ## called here with the same alpha and aggregate the fit will use --
        ## and K_default has no meaning for it.
        K_default <- NA_integer_
        n_flip_default <- .default_n_flip(NULL, dims[1L], dims[2L],
                                          alpha = alpha, aggregate = aggregate)
      }
    } else if (design == "panel") {
      if (is.null(time_v) && C != 3L)
        stop(paste0("design = \"panel\" needs a time dimension: pass 3 ",
                    "index columns or `time =`."), call. = FALSE)
      third_index_as_time()
      pc <- panel_cells()
      if (!pc$complete)
        stop(sprintf(paste0(
          "design = \"panel\" needs a complete balanced (row, col, time) ",
          "array, but only %d of %s = %d cells are observed. Use design = ",
          "\"panel_missing\" (or call mwperm_panel_missing() directly), ",
          "which restricts to the (row, col) pairs observed in every period ",
          "and permutes within fully observed blocks with the period held ",
          "fixed; or curate a complete balanced subset."),
          N, pc$dims_str, pc$expected), call. = FALSE)
      roles <- list(row = names(idx)[1L], col = names(idx)[2L],
                    time = names(time_v))
      K_default <- min(dims[1:2]) - 1L
      cells_obs <- N
      cells_exp <- pc$expected
      balance <- "complete"
    } else if (design == "panel_missing") {
      if (is.null(time_v) && C != 3L)
        stop(paste0("design = \"panel_missing\" needs a time dimension: ",
                    "pass 3 index columns or `time =`."), call. = FALSE)
      third_index_as_time()
      pc <- panel_cells()
      finish_panel_missing("forced via design =", complete = pc$complete,
                           expected = pc$expected, announce = FALSE)
    } else if (design == "threeway") {
      if (C != 3L)
        stop("design = \"threeway\" needs exactly 3 index dimensions.",
             call. = FALSE)
      cells3 <- .cell_code(cbind(dense[[1L]], dense[[2L]], dense[[3L]]))
      cells_obs <- N
      cells_exp <- prod(dims)
      if (anyDuplicated(cells3) > 0L || N != cells_exp)
        stop(sprintf(paste0(
          "design = \"threeway\" needs a complete balanced crossed array: %d ",
          "observations vs %s = %d expected cells%s. Curate a balanced subset ",
          "or use the 2-index dyadic/missing design."),
          N, paste(sprintf("%s=%d", names(dims), dims), collapse = " x "),
          cells_exp,
          if (anyDuplicated(cells3) > 0L) "; some cells repeat" else ""),
          call. = FALSE)
      roles <- list(id1 = names(idx)[1L], id2 = names(idx)[2L],
                    id3 = names(idx)[3L])
      K_default <- min(dims) - 1L
      balance <- "complete"
      ## Forcing threeway permutes EVERY index: if one looks
      ## like time -- by name, or strongly by value (temporal class, or
      ## regularly spaced with strictly fewest levels) -- say so, because a
      ## permuted time series over-rejects badly under serial dependence.
      tl <- union(which(vapply(names(idx), .timelike_name, logical(1))),
                  which(vapply(seq_len(3L), function(k)
                    .timelike_strong(idx[[k]], dims[k], min(dims[-k])),
                    logical(1))))
      if (length(tl))
        warns <- c(warns, paste0(
          "design = \"threeway\" permutes every index, but ",
          paste0("'", names(idx)[tl], "'", collapse = ", "),
          if (length(tl) > 1L) " look" else " looks",
          " time-like. If any of these is a time dimension the three-way ",
          "test is invalid under serial dependence (it can over-reject ",
          "badly): use design = \"panel\" with time = naming the time ",
          "index. Proceed only if the errors are exchangeable in every ",
          "index."))
    } else if (design == "layout") {
      if (C != 2L)
        stop("design = \"layout\" needs exactly 2 index dimensions (cells).",
             call. = FALSE)
      finish_layout("forced via design =")
      reason <- "forced via design ="
    } else if (design == "irregular" && !is.null(time_v)) {
      ## Section 6.4 with a PERIOD as the within-cell index is the incomplete
      ## panel with one observation per cell and period (the advisors'
      ## reading, 2026-09-22), and the random per-cell trim is exact only for
      ## exchangeable replicates: it keeps different periods in different
      ## cells and can over-reject under a period effect even with a
      ## cell-constant d. So a tagged `time =` runs mwperm_panel_missing(),
      ## which keeps the same L0 periods in every cell and holds the period
      ## fixed. (Until 0.4.2 `time` was dropped here without a word and the
      ## random trim ran on the order of appearance.)
      if (C != 2L)
        stop(paste0("design = \"irregular\" with `time =` needs exactly 2 ",
                    "index dimensions (the cells)."), call. = FALSE)
      if (!is.null(rep_v))
        stop(paste0("design = \"irregular\" takes one within-cell index: ",
                    "`time =` if it is a period (the fit then runs ",
                    "mwperm_panel_missing()), `rep =` if the repeats are ",
                    "exchangeable replicates -- not both."), call. = FALSE)
      pc <- panel_cells()
      finish_panel_missing(paste0("design = \"irregular\" with a `time =` ",
                                  "index: periods run as the incomplete ",
                                  "panel"),
                           complete = pc$complete, expected = pc$expected,
                           announce = FALSE)
      notes <- c(notes, paste0(
        "design = \"irregular\" was given a `time =` index, so this is ",
        "mwperm_panel_missing(): the random per-cell trim of ",
        "mwperm_irregular() is exact only for exchangeable replicates, and ",
        "for periods the same L0 periods must be kept in every cell (pass ",
        "`L0 =` for that; without it, only the pairs observed in every ",
        "period are kept)."))
    } else if (design == "irregular") {
      if (C != 2L)
        stop(paste0("design = \"irregular\" needs exactly 2 index ",
                    "dimensions (the cells); the within-cell index is the ",
                    "`rep =` role."), call. = FALSE)
      if (!dup2)
        stop(paste0("design = \"irregular\" needs repeated (row, col) ",
                    "cells: Section 6.4 reduces each cell to L0 ",
                    "observations, which needs more than one per cell. With ",
                    "one observation per cell use design = \"dyadic\" (or ",
                    "\"missing\" if the array is incomplete)."),
             call. = FALSE)
      finish_irregular("forced via design =")
      reason <- "forced via design ="
    }
  } else if (!is.null(rep_v)) {
    ## user declared within-cell replication
    if (C != 2L)
      stop("With `rep =`, pass exactly 2 index dimensions (the cells).",
           call. = FALSE)
    finish_layout("`rep =` declares within-cell replication")
  } else if (!is.null(time_v)) {
    ## user declared a time dimension -> panel; an incomplete array is the
    ## incomplete-panel design, exactly as on the three-untagged-index route
    ## below (completeness is structure, so this fork is silent)
    if (C != 2L)
      stop(paste0("With `time =`, pass exactly 2 index dimensions (the ",
                  "cross-section)."), call. = FALSE)
    pc <- panel_cells()
    if (pc$complete) {
      chosen <- "panel"
      roles <- list(row = names(idx)[1L], col = names(idx)[2L],
                    time = names(time_v))
      reason <- "`time =` declares a panel"
      K_default <- min(dims[1:2]) - 1L
      cells_obs <- N
      cells_exp <- pc$expected
      balance <- "complete"
    } else {
      finish_panel_missing("`time =` declares a panel; incomplete array",
                           complete = FALSE, expected = pc$expected,
                           announce = TRUE)
    }
  } else if (C == 2L) {
    if (dup2) {
      ## ASSUMPTION FORK: layout by default, suppressed-panel warning
      finish_layout(
        "repeated (row, col) cells read as within-cell replication",
        paste0("Repeated cells were treated as exchangeable within-cell ",
               "replication (layout design). If the repeats are really time ",
               "periods, that assumption is false and the layout test is ",
               "invalid -- pass the time variable via `time =` to run the ",
               "panel test instead."))
    } else {
      cells_obs <- N
      cells_exp <- prod(dims[1:2])
      roles <- list(row = names(idx)[1L], col = names(idx)[2L])
      K_default <- min(dims[1:2]) - 1L
      if (N == cells_exp) {
        chosen <- "dyadic"
        balance <- "complete"
        reason <- "2 indices, one observation per cell, complete array"
        ## Heteroskedasticity leaves no trace in the clustering structure, so
        ## the sign-flip test can never be detected -- only offered.
        alternatives <- c(alternatives, paste0(
          "design = \"dyadic_het\" runs the sign-flip test instead: valid ",
          "under arbitrary heteroskedasticity for errors that are ",
          "independent across cells and symmetric about zero, but NOT under ",
          "additive cluster effects eta_i + xi_j (see ?mwperm_dyadic_het)"))
      } else {
        chosen <- "missing"
        balance <- sprintf("incomplete (%d of %d cells)", N, cells_exp)
        reason <- sprintf("2 indices, %d of %d cells observed", N, cells_exp)
        K_default <- NA_integer_          # depends on the biclique blocks
        notes <- c(notes, paste0(
          "The permutation-group order under missingness is set by the fully ",
          "observed biclique blocks; see find_bicliques() for the ",
          "achievable K."))
        ## The sign-flip test needs no complete array (nothing is discarded),
        ## which makes it worth naming here; its assumption is the same as
        ## on a complete array.
        alternatives <- c(alternatives, paste0(
          "design = \"dyadic_het\" runs the sign-flip test on every observed ",
          "cell, with no biclique search and nothing discarded: valid under ",
          "arbitrary heteroskedasticity for errors that are independent ",
          "across cells and symmetric about zero, but NOT under additive ",
          "cluster effects eta_i + xi_j (see ?mwperm_dyadic_het)"))
      }
    }
  } else {
    ## C == 3, no tags. A REPEATED cell is still fatal -- with two
    ## observations in one (i, j, t) neither the complete-array designs nor the
    ## blockwise one has a well-defined map. An INCOMPLETE array is not: since
    ## 0.4.0 it routes to mwperm_panel_missing(), which masks to the pairs
    ## observed in every period and permutes blockwise. The time role is
    ## decided first because that choice does not depend on completeness, and
    ## the incomplete design needs it too.
    cells3 <- .cell_code(cbind(dense[[1L]], dense[[2L]], dense[[3L]]))
    cells_obs <- N
    cells_exp <- prod(dims)
    if (anyDuplicated(cells3) > 0L)
      stop(sprintf(paste0(
        "3 index dimensions but some (%s) cell repeats, so no crossed design ",
        "applies: %d observations for %d cells. If the repeats are ",
        "replication within a cell, pass that index as `rep =`; otherwise ",
        "aggregate to one observation per cell."),
        paste(names(dims), collapse = ", "), N, cells_exp), call. = FALSE)
    ## which index is time? tagged > name > class/values > ambiguous
    name_hit <- which(vapply(names(idx), .timelike_name, logical(1)))
    val_hit <- which(vapply(seq_len(3L), function(k)
      .timelike_values(idx[[k]], dims[k], max(dims[-k])), logical(1)))
    strong_hit <- which(vapply(seq_len(3L), function(k)
      .timelike_strong(idx[[k]], dims[k], min(dims[-k])), logical(1)))
    t_k <- NULL
    t_why <- NULL
    if (length(name_hit) == 1L) {
      t_k <- name_hit
      t_why <- sprintf("'%s' identified as time by name", names(idx)[t_k])
      ## A name match alone is weak evidence: warn unless the
      ## values single this index out. Assigning the time role to a mere
      ## cluster sends the TRUE time dimension into the permuted pair, and a
      ## permuted time series over-rejects badly.
      if (!(length(strong_hit) == 1L && strong_hit == t_k)) {
        others <- setdiff(strong_hit, t_k)
        warns <- c(warns, paste0(
          "The time role went to '", names(idx)[t_k], "' on the strength of ",
          "its name alone",
          if (length(others))
            paste0(" -- and ", paste0("'", names(idx)[others], "'",
                                      collapse = ", "),
                   if (length(others) > 1L) " look" else " looks",
                   " at least as time-like by value")
          else " (its values do not single it out)",
          ". If the true time dimension is a different index, permuting it ",
          "makes the test invalid (it can over-reject badly). Tag the ",
          "correct index explicitly via time = ; force design = \"threeway\" ",
          "only if all three dimensions are exchangeable."))
      }
    } else if (length(val_hit) == 1L) {
      t_k <- val_hit
      t_why <- sprintf(paste0("'%s' identified as time by its values ",
                              "(temporal/regularly spaced)"),
                       names(idx)[t_k])
    }
    if (is.null(t_k)) {
      ## ASSUMPTION FORK: ambiguous -> panel default on the third index
      t_k <- 3L
      t_why <- sprintf("ambiguous third dimension; '%s' held fixed by default",
                       names(idx)[3L])
      warns <- c(warns, paste0(
        "No index is clearly the time dimension, so the PANEL design was ",
        "chosen with '", names(idx)[3L], "' held fixed. That default is ",
        "protective ONLY if the permuted pair ('", names(idx)[1L], "', '",
        names(idx)[2L], "') is exchangeable -- if one of THOSE is the true ",
        "time dimension the test is invalid, so tag the correct index via ",
        "time = . If all three dimensions really are exchangeable, force ",
        "design = \"threeway\" to recover power."))
    }
    time_v <- idx[t_k]
    idx <- idx[-t_k]
    dims_cs <- dims[-t_k]
    roles <- list(row = names(idx)[1L], col = names(idx)[2L],
                  time = names(time_v))
    reason <- t_why
    if (N == cells_exp) {
      chosen <- "panel"
      K_default <- min(dims_cs) - 1L
      balance <- "complete"
    } else {
      finish_panel_missing(paste0(t_why, "; incomplete array"),
                           complete = FALSE, expected = cells_exp,
                           announce = TRUE)
    }
    dims <- c(dims_cs, dims[t_k])
    names(dims) <- c(names(idx), names(time_v))
  }

  K_default <- if (is.na(K_default)) NA_integer_ else min(K_default, 199L)
  ## The verdict is about the p-value the FIT will report, so it uses the same
  ## floor .ipt_engine() gates on: 1/(K+1) under "median", 2/(K+1) under
  ## "median2" (which reports min(1, 2 x median)), compared against the alpha
  ## the fit will use -- not a hard-coded 0.05. For the sign-flip design the
  ## group order is 2^(n_flip - 1) at the default n_flip, and `levels_needed`
  ## is then the n_flip a set needs rather than a cluster count.
  agg_mult <- if (identical(aggregate, "median2")) 2L else 1L
  need_lvl <- as.integer(ceiling(agg_mult / alpha))
  if (!is.null(n_flip_default)) {
    p_floor <- min(1, agg_mult / 2^(n_flip_default - 1L))
    need_lvl <- 1L + as.integer(ceiling(log2(need_lvl)))
    res_ok <- p_floor <= alpha
  } else {
    p_floor <- if (is.na(K_default)) NA_real_ else
      min(1, agg_mult / (K_default + 1L))
    res_ok <- if (is.na(K_default)) NA else p_floor <= alpha
  }

  ## the exact downstream call
  fn <- c(dyadic = "mwperm_dyadic", panel = "mwperm_panel",
          panel_missing = "mwperm_panel_missing",
          threeway = "mwperm_threeway", layout = "mwperm_layout",
          missing = "mwperm_missing",
          irregular = "mwperm_irregular",
          dyadic_het = "mwperm_dyadic_het")[[chosen]]
  args <- switch(chosen,
    dyadic = sprintf("row = %s, col = %s", roles$row, roles$col),
    dyadic_het = sprintf("row = %s, col = %s", roles$row, roles$col),
    missing = sprintf("row = %s, col = %s, min_block = ...", roles$row,
                      roles$col),
    panel = sprintf("row = %s, col = %s, time = %s, time_fe = TRUE",
                    roles$row, roles$col, roles$time),
    panel_missing = sprintf(paste0("row = %s, col = %s, time = %s, ",
                                   "L0 = NULL, min_block = ..., ",
                                   "time_fe = TRUE"),
                            roles$row, roles$col, roles$time),
    threeway = sprintf("id1 = %s, id2 = %s, id3 = %s",
                       roles$id1, roles$id2, roles$id3),
    layout = sprintf("row = %s, col = %s%s", roles$row, roles$col,
                     if (identical(roles$rep, "(within-cell order)")) ""
                     else sprintf(", rep = %s", roles$rep)),
    irregular = sprintf("row = %s, col = %s%s, L0 = ...", roles$row, roles$col,
                        if (identical(roles$rep, "(within-cell order)")) ""
                        else sprintf(", rep = %s", roles$rep)))
  call_str <- sprintf("%s(y, d, x, %s)", fn, args)

  structure(list(
    design = chosen, roles = roles, dims = dims, n_obs = N,
    cells = c(observed = cells_obs, expected = cells_exp),
    balance = balance, K_default = K_default,
    alpha = alpha, aggregate = aggregate, p_floor = p_floor,
    levels_needed = need_lvl, resolution_ok = res_ok,
    n_flip_default = n_flip_default,
    call_str = call_str, reason = reason,
    warnings = warns, notes = notes, alternatives = alternatives,
    time = if (!is.null(time_v)) time_v[[1L]] else NULL,
    rep = if (!is.null(rep_v)) rep_v[[1L]] else NULL,
    index = idx
  ), class = "mwperm_design")
}

#' @rdname mwperm_check
#' @param x An object of class `"mwperm_design"` (print method).
#' @param ... Ignored.
#' @export
print.mwperm_design <- function(x, ...) {
  cat("\nmwperm design diagnosis\n")
  cat(strrep("-", 30), "\n", sep = "")
  cat("Detected design : ", x$design, " (", x$reason, ")\n", sep = "")
  role_txt <- paste(sprintf("%s = %s", names(x$roles), unlist(x$roles)),
                    collapse = ", ")
  cat("Roles           : ", role_txt, "\n", sep = "")
  cat("Dimensions      : ",
      paste(sprintf("%s (%d)", names(x$dims), x$dims), collapse = " x "),
      sprintf(" | %d observations\n", x$n_obs), sep = "")
  cat("Balance         : ", x$balance, "\n", sep = "")
  is_flip <- !is.null(x$n_flip_default)
  if (is.na(x$K_default) && !is_flip) {
    cat("Resolution      : set by the biclique blocks, so not known until\n",
        "                  they are found (see find_bicliques)\n", sep = "")
  } else {
    ## Report the resolution as both the fraction and its decimal: the
    ## fraction shows where it comes from, the decimal is what gets compared
    ## against alpha. "Too coarse" means no 95% set can exclude anything --
    ## the p-value stays exact either way, so say what is and is not lost.
    ## The sign-flip design's order is 2^(n_flip - 1), and its remedy is a
    ## larger n_flip, not more clusters.
    order <- if (is_flip) 2^(x$n_flip_default - 1L) else x$K_default + 1L
    if (is_flip)
      cat(sprintf(paste0("Resolution      : default n_flip = %d, so p-values ",
                         "are multiples of 1/2^%d = 1/%d = %s\n"),
                  x$n_flip_default, x$n_flip_default - 1L, order,
                  .fmt_p(1 / order)))
    else
      cat(sprintf(paste0("Resolution      : default K = %d, so p-values are ",
                         "multiples of 1/%d = %s\n"),
                  x$K_default, order, .fmt_p(1 / order)))
    ## Objects from before `alpha`/`aggregate` were fields carry neither; read
    ## the 0.05 / "median" they were computed under.
    alpha <- if (is.null(x$alpha)) 0.05 else x$alpha
    agg <- if (is.null(x$aggregate)) "median" else x$aggregate
    need <- if (is.null(x$levels_needed)) 20L else x$levels_needed
    lvl <- sprintf("%.0f%%", 100 * (1 - alpha))
    art <- if (substr(lvl, 1L, 1L) == "8") "an" else "a"   # "an 80%" set
    verdict <- if (isTRUE(x$resolution_ok))
      sprintf("-> fine enough for %s %s confidence set at alpha = %s",
              art, lvl, format(alpha))
    else sprintf(paste0("-> TOO COARSE for %s %s confidence set at alpha = ",
                        "%s (p cannot reach %s). The p-value is still ",
                        "exact; %s %s set needs %s."),
                 art, lvl, format(alpha), format(alpha), art, lvl,
                 if (is_flip) sprintf("n_flip >= %d", need)
                 else sprintf(">= %d levels in the smallest permuted dimension",
                              need))
    if (identical(agg, "median2"))
      verdict <- paste0(verdict, sprintf(paste0(
        " Under aggregate = \"median2\" the reported p-value is min(1, 2 x ",
        "median), so its floor is 2/%d = %s."),
        order, .fmt_p(2 / order)))
    ## Wrapped at the full console width (not the 0.9 the notes use) so the
    ## verdict's first line -- the one README.md shows -- stays whole.
    cat(strwrap(verdict, initial = "                  ",
                prefix = "                     ",
                width = getOption("width", 80)), sep = "\n")
  }
  cat("Would run       : ", x$call_str, "\n", sep = "")
  for (w in x$warnings)
    cat(strwrap(w, initial = "  ! ", prefix = "    ",
                width = 0.9 * getOption("width", 80)), sep = "\n")
  for (nt in x$notes)
    cat(strwrap(nt, initial = "  - ", prefix = "    ",
                width = 0.9 * getOption("width", 80)), sep = "\n")
  ## Designs the data cannot select (the assumption is invisible in the
  ## structure), offered in one line each. Objects fitted before this field
  ## existed carry NULL, which the loop treats as none.
  for (alt in x$alternatives)
    cat(strwrap(alt, initial = "  ? ", prefix = "    ",
                width = 0.9 * getOption("width", 80)), sep = "\n")
  cat("\n")
  invisible(x)
}

## ---- mwperm: one-call dispatch ----------------------------------------------

#' One-call invariant permutation test with automatic design detection
#'
#' Detects the clustering design of the data via [mwperm_check()] and
#' dispatches to the matching test -- [mwperm_dyadic()], [mwperm_panel()],
#' [mwperm_panel_missing()], [mwperm_threeway()], [mwperm_layout()] or
#' [mwperm_missing()]; [mwperm_irregular()] is never chosen automatically and
#' needs `design = "irregular"` with `L0` (with a `time =` index that runs
#' [mwperm_panel_missing()] instead), and [mwperm_dyadic_het()] is never
#' chosen automatically and needs `design = "dyadic_het"` -- forwarding all
#' arguments unchanged. A thin convenience layer: the returned object is
#' exactly what the underlying function returns (plus a record of what was
#' detected), and calling the specific function directly with the same seed
#' gives identical results.
#'
#' See [mwperm_check()] for the detection rules, in particular the two
#' assumption-dependent forks (panel-vs-threeway and
#' layout-vs-suppressed-panel) that are announced rather than silently
#' resolved. Structural forks (complete vs incomplete arrays, replicated
#' cells) are resolved silently.
#'
#' @param y,d,x Outcome, covariate(s) of interest, and optional nuisance
#'   covariates, as in [mwperm_dyadic()]. With `data` given, each may also be
#'   a character (vector of) column name(s) resolved against it.
#' @param index The clustering dimensions (2 or 3): a data frame, named list
#'   of vectors, or character vector of column names in `data`.
#' @param data Optional data frame; column names in `y`, `d`, `x`, `index`,
#'   `time`, `rep` are resolved against it.
#' @param time,rep Optional explicit role tags (vector or column name); see
#'   [mwperm_check()].
#' @param design Force a design instead of auto-detecting (the structure is
#'   still validated against it). `"dyadic_het"` runs [mwperm_dyadic_het()],
#'   the sign-flip test that tolerates heteroskedasticity at the price of
#'   assuming errors independent across cells and symmetric (no additive
#'   cluster effects); it is *never* chosen automatically, because
#'   heteroskedasticity is invisible in the clustering structure. It runs
#'   on incomplete arrays too (nothing is discarded).
#' @param K Number of non-identity permutations; the default and the
#'   admissible range depend on the dispatched design -- see the dispatched
#'   function. Does not apply to `design = "dyadic_het"`, whose group is
#'   sized by `n_flip`; supplying it there warns and ignores it.
#' @param n_flip Passed to [mwperm_dyadic_het()] (`design = "dyadic_het"`
#'   only: the number of flip groups, group order `2^(n_flip - 1)`; by
#'   default the smallest number whose p-value floor `1 / 2^(n_flip - 1)`
#'   (doubled under `aggregate = "median2"`) is at most `alpha`, i.e. 6 at
#'   `alpha = 0.05`. Supplying it for another design warns and ignores it).
#' @param time_fe Passed to [mwperm_panel()] or [mwperm_panel_missing()]
#'   (panel designs only; supplying it for another design warns and ignores
#'   it).
#' @param L0 Passed to [mwperm_layout()], [mwperm_irregular()] or
#'   [mwperm_panel_missing()] (layout, irregular and incomplete-panel designs
#'   only; required for `design = "irregular"`, optional for the other two:
#'   for an incomplete panel it keeps the L0 periods jointly observed by the
#'   most pairs, the same ones in every cell).
#' @param trim Removed in 0.4.2 and refused when supplied: the level-aligned
#'   cut that was `mwperm_irregular(trim = "levels")` is now
#'   `mwperm_panel_missing(time = <rep>, L0 = <L0>)` (reached here with
#'   `time =` and `L0 =`), and [mwperm_irregular()] implements only the
#'   paper's random per-cell trim.
#' @param min_block,block_method,permute `min_block` and `block_method` are
#'   passed to [mwperm_missing()], [mwperm_panel_missing()] or
#'   [mwperm_irregular()]; `permute` to [mwperm_missing()] only. Supplying any
#'   of them for another design warns and ignores it.
#' @param n_reps Number of repetitions whose p-values are aggregated (see
#'   `aggregate`). `NULL` (the default) means the dispatched function's own
#'   default: 500 for [mwperm_irregular()], whose random per-cell trim is
#'   redrawn in every repetition, and 10 for every other design. A value
#'   given here is forwarded as it is.
#' @param verbose If `TRUE` (default) print one line stating the detected
#'   design and the dispatched call.
#' @inheritParams mwperm_dyadic
#'
#' @param aggregate How the `n_reps` per-repetition p-values are combined into
#'   the reported p-value, and into the confidence set that inverts it.
#'   `"median"` (the default) is the median, as recommended in Remark 1 of
#'   Guo, Toulis and Wang (2026); `"median2"` is `min(1, 2 * median)`.
#'
#' The choice decides what "exact" covers. Theorem 1 gives finite-sample
#' validity for a single random permutation group, so at `n_reps = 1` the
#' p-value is exact as stated. The median of several dependent randomised
#' p-values is a de-randomisation heuristic: endorsed by Remark 1 and well
#' behaved in practice, but not itself guaranteed valid at level `alpha`.
#' Twice the median is guaranteed, under arbitrary dependence across
#' repetitions (Ruschendorf 1982; Vovk and Wang 2020).
#'
#' So use `"median2"` when the guarantee must hold as stated with `n_reps >
#' 1`. It is conservative: it never rejects where `"median"` would not, and
#' its confidence set is never narrower. The default is unchanged, so existing
#' numbers stand.
#'
#' The cost is resolution. `"median2"` reports `min(1, 2 * median)`, so its
#' smallest attainable p-value is `2/(K+1)`, not `1/(K+1)`, and rejecting at
#' level `alpha` needs `K + 1 >= 2/alpha` -- at `alpha = 0.05` that is 40
#' levels in the smallest permuted dimension, twice what `"median"` needs.
#' Below that the p-value is still exact but cannot reach `alpha`, and the fit
#' says so in a note.
#' @return The `"mwperm"` object of the dispatched test, with an extra `auto`
#'   field recording the detection (design, roles, reason); the detection
#'   notices are prepended to the object's `note` field and shown by
#'   [print.mwperm()], and any assumption-fork or weak-evidence notice is
#'   additionally raised as a `warning` at fit time. Field provenance (see
#'   [mwperm_dyadic()] for the full account): `estimate`/`se_naive` are the
#'   OLS estimate and naive SE, `conf_int` (or `conf_region`/`conf_box` for
#'   several coefficients) the IPT inverted-test confidence set, and `pvalue`
#'   the IPT permutation p-value.
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data. arXiv:2601.08610.
#' @seealso [mwperm_check()] for the diagnosis without any computation;
#'   [mwperm_dyadic()], [mwperm_panel()], [mwperm_threeway()],
#'   [mwperm_layout()], [mwperm_irregular()], [mwperm_missing()],
#'   [mwperm_dyadic_het()] for the underlying tests.
#' @examples
#' data(trade_dyadic)
#' fit <- mwperm(y = "log_trade", d = "log_dist",
#'               x = c("log_gdp_i", "log_gdp_j"),
#'               index = c("importer", "exporter"),
#'               data = trade_dyadic, seed = 1)
#' fit
#' @export
mwperm <- function(y, d, x = NULL, index, data = NULL, time = NULL, rep = NULL,
                   design = c("auto", "dyadic", "threeway", "panel",
                              "panel_missing", "layout", "missing",
                              "irregular", "dyadic_het"),
                   K = NULL, alpha = 0.05, beta_null = 0, conf_int = TRUE,
                   n_reps = NULL, seed = NULL, grid = NULL, n_cores = 1L,
                   time_fe = TRUE, L0 = NULL, trim = NULL,
                   min_block = 3L, block_method = c("greedy", "exact"),
                   permute = c("both", "rows", "cols"), n_flip = NULL,
                   aggregate = c("median", "median2"), verbose = TRUE) {
  design <- match.arg(design)
  aggregate <- match.arg(aggregate)
  ## `trim` is kept as a slot for one release so that a 0.4.1 call fails with
  ## the replacement named rather than with R's "unused argument".
  if (!is.null(trim))
    stop(paste0("`trim` was removed in 0.4.2: for periods with a common ",
                "effect use mwperm_panel_missing(time = <rep>, L0 = <L0>) ",
                "(through mwperm(), pass the period as `time =` together ",
                "with `L0 =`); mwperm_irregular() now implements only the ",
                "paper's random per-cell trim."), call. = FALSE)
  cl <- match.call()
  ## capture the caller's expression for d BEFORE evaluation: the front ends
  ## label coefficients by deparse(substitute(d)), which through do.call would
  ## deparse the VALUES; we forward d as a named-column matrix instead.
  ## A multi-element deparse means exactly that -- the argument arrived as a
  ## value, not an expression -- so collapsing it would make the whole data
  ## vector the coefficient label. Use a generic name instead (see
  ## .coef_names, which applies the same rule to the direct front ends).
  d_expr <- deparse(substitute(d))
  d_expr <- if (length(d_expr) == 1L) d_expr else "d"

  ## resolve y / d / x against `data`: accept bare vectors, column-name
  ## strings, or character vectors of column names (for x)
  resolve_col <- function(v, what, allow_multi = FALSE) {
    if (is.null(v)) return(NULL)
    if (is.character(v) && !is.null(data) && all(v %in% names(data))) {
      if (length(v) == 1L && !allow_multi) {
        out <- data[[v]]
        attr(out, "mwperm_name") <- v
        return(out)
      }
      if (allow_multi) {
        out <- as.matrix(data[v])
        return(out)
      }
    }
    v
  }
  y <- resolve_col(y, "y")
  d <- resolve_col(d, "d", allow_multi = FALSE)
  if (is.character(d) && !is.null(data) && all(d %in% names(data)))
    d <- as.matrix(data[d])                      # multi-column d by names
  x <- resolve_col(x, "x", allow_multi = TRUE)
  for (nm in c("y", "d", "x")) {                 # names that failed to resolve
    v <- get(nm)
    if (is.character(v))
      stop(sprintf(paste0("`%s` is a character vector; to use column ",
                          "names, pass ",
                          "`data` containing column(s) %s."),
                   nm, paste0("'", v, "'", collapse = ", ")), call. = FALSE)
  }

  ## carry a readable coefficient label through the dispatch (colnames take
  ## precedence in .coef_names, so naming here fully avoids the deparse path)
  d_lab <- attr(d, "mwperm_name") %||% d_expr
  attr(d, "mwperm_name") <- NULL
  if (is.null(dim(d))) {
    d <- matrix(as.numeric(d), ncol = 1L, dimnames = list(NULL, d_lab))
  } else if (is.null(colnames(d))) {
    colnames(d) <- if (ncol(d) == 1L) d_lab else paste0(d_lab, seq_len(ncol(d)))
  }
  attr(y, "mwperm_name") <- NULL

  chk <- mwperm_check(index = index, y = y, d = d, data = data,
                      time = time, rep = rep, design = design,
                      alpha = alpha, aggregate = aggregate)

  ## design-specific arguments must not be silently accepted for the wrong
  ## design (warn and ignore, mirroring the printed diagnosis)
  supplied <- names(cl)
  check_arg <- function(arg, ok_design) {
    nd <- length(ok_design)
    listed <- if (nd <= 2L) paste(ok_design, collapse = " and ") else
      paste0(paste(ok_design[-nd], collapse = ", "), " and ", ok_design[nd])
    if (arg %in% supplied && !chk$design %in% ok_design)
      warning(sprintf(paste0("`%s` applies to the %s design%s only; it was ",
                             "ignored for '%s'."),
                      arg, listed, if (nd > 1L) "s" else "",
                      chk$design), call. = FALSE)
  }
  check_arg("time_fe", c("panel", "panel_missing"))
  check_arg("L0", c("layout", "irregular", "panel_missing"))
  check_arg("min_block", c("missing", "irregular", "panel_missing"))
  check_arg("block_method", c("missing", "irregular", "panel_missing"))
  check_arg("permute", "missing")
  check_arg("n_flip", "dyadic_het")
  ## The sign-flip group is sized by n_flip, not K (its order is 2^(n_flip -
  ## 1)); a K passed to it would be silently meaningless, so say so.
  if ("K" %in% supplied && chk$design == "dyadic_het")
    warning(paste0("`K` sizes a permutation group and does not apply to the ",
                   "sign-flip design; it was ignored for 'dyadic_het'. Use ",
                   "`n_flip` (group order 2^(n_flip - 1))."), call. = FALSE)
  if (chk$design == "irregular" && is.null(L0))
    stop(paste0("The Section 6.4 (irregular) design requires `L0 =`, the ",
                "number of observations retained per cell, which sets the ",
                "mask M_ij = 1{ell_ij >= L0}. See ?mwperm_irregular."),
         call. = FALSE)

  ## Assumption-fork and weak-evidence detection notices are REAL warnings at
  ## fit time -- they flag branches that can be anti-conservative if the
  ## guessed role is wrong, and must be impossible to miss.
  ## They also stay on the returned object's `note`, as before.
  for (w in chk$warnings) warning(w, call. = FALSE)

  if (isTRUE(verbose)) {
    ## Two lines. The single-line form ran to ~140 characters and wrapped
    ## wherever the terminal happened to end, splitting the dispatched call
    ## mid-argument; this is also the shape README.md documents.
    message(sprintf("Detected design: %s (%s)", chk$design, chk$reason),
            "\n  -> running ", chk$call_str)
  }

  ## n_reps is forwarded only when given: NULL means the dispatched front
  ## end's own default, which is 500 for mwperm_irregular() (its random
  ## per-cell trim is redrawn every repetition) and 10 for every other design.
  common <- c(list(y = y, d = d, x = x, K = K, alpha = alpha,
                   beta_null = beta_null, conf_int = conf_int),
              if (!is.null(n_reps)) list(n_reps = n_reps),
              list(seed = seed, grid = grid, aggregate = aggregate,
                   n_cores = n_cores))
  ix <- chk$index
  res <- switch(chk$design,
    dyadic = do.call(mwperm_dyadic,
                     c(common, list(row = ix[[1L]], col = ix[[2L]]))),
    dyadic_het = do.call(mwperm_dyadic_het,
                         c(common[names(common) != "K"],
                           list(row = ix[[1L]], col = ix[[2L]],
                                n_flip = n_flip))),
    missing = do.call(mwperm_missing,
                      c(common, list(row = ix[[1L]], col = ix[[2L]],
                                     min_block = min_block,
                                     block_method = block_method,
                                     permute = permute))),
    panel = do.call(mwperm_panel,
                    c(common, list(row = ix[[1L]], col = ix[[2L]],
                                   time = chk$time, time_fe = time_fe))),
    panel_missing = do.call(mwperm_panel_missing,
                            c(common, list(row = ix[[1L]], col = ix[[2L]],
                                           time = chk$time, L0 = L0,
                                           min_block = min_block,
                                           block_method = block_method,
                                           time_fe = time_fe))),
    threeway = do.call(mwperm_threeway,
                       c(common, list(id1 = ix[[1L]], id2 = ix[[2L]],
                                      id3 = ix[[3L]]))),
    layout = do.call(mwperm_layout,
                     c(common, list(row = ix[[1L]], col = ix[[2L]],
                                    rep = chk$rep, L0 = L0))),
    irregular = do.call(mwperm_irregular,
                        c(common, list(row = ix[[1L]], col = ix[[2L]],
                                       rep = chk$rep, L0 = L0,
                                       min_block = min_block,
                                       block_method = block_method))))

  ## make the automatic choice transparent on the returned object
  res$auto <- list(design = chk$design, reason = chk$reason, roles = chk$roles)
  res$note <- c(chk$warnings, chk$notes, res$note)
  res$call <- cl
  res
}
