## Unified entry point: automatic design detection (mwperm_check) and
## dispatch (mwperm). A thin, additive layer over the five mwperm_* front
## ends -- it changes nothing about how any test is computed.
##
## Detection policy (see REVIEW.md 8 and the decision table in
## docs/REVIEW_TASK.md): forks that are visible in the data's *structure*
## (index count, repeated cells, completeness) are resolved silently; forks
## that hinge on an EXCHANGEABILITY ASSUMPTION the data cannot reveal are
## announced, default to the choice that stays valid under the widest set of
## data-generating processes, and carry override instructions:
##   * 3 complete crossed indices: panel vs threeway is such a fork. Running
##     threeway on a panel (time-autocorrelated errors) is INVALID (size
##     distortion; empirically 0.88 rejection at alpha = .2 on a trending
##     panel null -- REVIEW.md 6); running panel on genuinely three-way
##     exchangeable data is merely less powerful (InvA implies InvB). Hence
##     panel is the default; threeway only on explicit design = "threeway".
##   * 2 indices with repeated cells: layout vs a panel whose time index was
##     not passed. Layout assumes the within-cell replicates are exchangeable;
##     if they are a time series that is false. Default layout + warning.

## ---- small helpers ----------------------------------------------------------

#' Resolve the `index` argument to a named list of equal-length vectors.
#' @keywords internal
#' @noRd
.resolve_index <- function(index, data) {
  if (is.character(index) && is.null(dim(index))) {
    ## character vector of column names against `data`
    if (is.null(data))
      stop("`index` is a character vector of column names but `data` was not supplied.",
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

## ---- mwperm_check -----------------------------------------------------------

#' Diagnose a multi-way clustered dataset and choose the appropriate design
#'
#' Inspects the clustering structure of a dataset -- number of index
#' dimensions, repeated cells, completeness/balance -- and reports which
#' \code{mwperm_*} test applies, without running any permutations or fitting
#' anything. \code{\link{mwperm}} uses it for automatic dispatch; call it
#' directly to see the diagnosis.
#'
#' Structural forks (dyadic vs missing vs layout-by-replication) are resolved
#' silently from the data. Two forks depend on an \emph{exchangeability
#' assumption the data cannot reveal} and are therefore announced with
#' override instructions, defaulting to the choice that remains valid under
#' the widest set of error processes:
#' \itemize{
#'   \item \strong{panel vs three-way} (complete balanced 3-index arrays):
#'     running \code{\link{mwperm_threeway}} on a panel whose errors are
#'     dependent over time is \emph{invalid}, while running
#'     \code{\link{mwperm_panel}} on genuinely three-way exchangeable data is
#'     valid, only less powerful. The default is therefore \code{panel}; a
#'     time-like third index (by name, temporal class, or regular numeric
#'     spacing) is used as the time dimension silently, an ambiguous one with
#'     a notice. Force \code{design = "threeway"} only when all three
#'     dimensions are genuinely exchangeable.
#'   \item \strong{layout vs suppressed panel} (2 indices with repeated
#'     cells): repeats are treated as within-cell replication
#'     (\code{\link{mwperm_layout}}), which assumes the replicates are
#'     exchangeable within cells -- if they are really a time series, pass the
#'     time variable via \code{time =} to get the panel test instead. A
#'     notice is attached.
#' }
#'
#' @param index The clustering dimensions (2 or 3): a data frame, a named
#'   list of vectors, or a character vector of column names resolved against
#'   \code{data}.
#' @param y,d Optional outcome and covariate(s) of interest; only used for
#'   extra diagnostics (e.g. the layout no-power warning when \code{d} is
#'   constant within every cell), never for fitting.
#' @param data Optional data frame against which character \code{index},
#'   \code{time} and \code{rep} entries are resolved.
#' @param time Optional explicit time dimension: a vector, or the name of a
#'   column of \code{data} (or of one of the \code{index} columns). Forces
#'   the panel interpretation of that dimension.
#' @param rep Optional explicit replication identifier (vector or column
#'   name): declares within-cell replication and forces the layout design.
#' @param design Force a design instead of auto-detecting (the structure is
#'   still validated against it).
#'
#' @return An object of class \code{"mwperm_design"}: a list with fields
#'   \code{design} (the chosen design), \code{roles} (which index plays
#'   row/col/id1..3/time/rep), \code{dims} (levels per dimension),
#'   \code{n_obs}, \code{cells} (observed/expected), \code{balance},
#'   \code{K_default} and \code{resolution_ok} (whether a 95\% confidence set
#'   is attainable), \code{call_str} (the downstream call), \code{reason}
#'   (one-line explanation), and \code{warnings}/\code{notes} (the
#'   assumption-fork notices etc.). Its \code{print} method lays this out as
#'   a short human diagnosis.
#'
#' @seealso \code{\link{mwperm}} for one-call dispatch.
#' @examples
#' data(trade_dyadic)
#' mwperm_check(index = c("importer", "exporter"), data = trade_dyadic)
#' data(trade_panel)
#' mwperm_check(index = c("importer", "exporter", "year"), data = trade_panel)
#' @export
mwperm_check <- function(index, y = NULL, d = NULL, data = NULL,
                         time = NULL, rep = NULL,
                         design = c("auto", "dyadic", "threeway", "panel",
                                    "layout", "missing")) {
  design <- match.arg(design)
  idx <- .resolve_index(index, data)
  N <- length(idx[[1L]])
  if (any(vapply(idx, length, integer(1)) != N))
    stop("All index dimensions must have the same length.", call. = FALSE)

  notes <- character(0); warns <- character(0)

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
      stop(sprintf("`%s = \"%s\"` matches neither an index column nor a column of `data`.",
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
      "Index dimension(s) %s have a single level and were dropped (no clustering).",
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
  dense <- lapply(idx, .dense_id)
  cells2 <- .cell_code(cbind(dense[[1L]], dense[[2L]]))
  dup2 <- anyDuplicated(cells2) > 0L

  ## defaults filled per design below
  roles <- NULL; chosen <- NULL; reason <- NULL
  K_default <- NA_integer_; balance <- NA_character_
  cells_obs <- NA_integer_; cells_exp <- NA_integer_

  finish_layout <- function(why, warn_txt = NULL) {
    cell <- .dense_id(interaction(dense[[1L]], dense[[2L]], drop = TRUE))
    sizes <- tabulate(cell)
    chosen <<- "layout"
    roles <<- list(row = names(idx)[1L], col = names(idx)[2L],
                   rep = if (!is.null(rep_v)) names(rep_v) else
                     "(within-cell order)")
    reason <<- why
    K_default <<- min(sizes) - 1L
    balance <<- sprintf("replicated (%d cells, %d-%d replicates)",
                        length(sizes), min(sizes), max(sizes))
    cells_obs <<- length(sizes); cells_exp <<- prod(dims[1:2])
    if (!is.null(warn_txt)) warns <<- c(warns, warn_txt)
    ## no-power diagnostic when d is available
    if (!is.null(d)) {
      D <- as.matrix(d)
      wv <- tapply(seq_len(N), cell, function(ii)
        any(apply(D[ii, , drop = FALSE], 2L, function(v) diff(range(v)) > 0)))
      if (!any(unlist(wv)))
        warns <<- c(warns, paste0(
          "`d` is constant within every (row, col) cell: the within-cell layout ",
          "test would have NO power. If the repeats are time periods, pass ",
          "`time =`; if `d` is a dyad-level covariate, this design cannot test it."))
    }
  }

  ## Completeness gate shared by every panel path (auto, tagged, forced):
  ## mwperm_panel hard-requires a complete balanced (row, col, time) array,
  ## so fail here with the actionable message rather than deep in the engine.
  require_complete_panel <- function() {
    nT <- length(unique(time_v[[1L]]))
    tri <- .cell_code(cbind(dense[[1L]], dense[[2L]], .dense_id(time_v[[1L]])))
    exp3 <- prod(dims[1:2]) * nT
    if (anyDuplicated(tri) > 0L || N != exp3)
      stop(sprintf(paste0(
        "A panel needs a complete balanced (row, col, time) array: %d ",
        "observations vs %s=%d x %s=%d x %s=%d = %d expected cells%s. ",
        "Options: (a) curate a complete balanced subset; (b) if the repeats ",
        "are exchangeable replication rather than time, pass them as `rep =` ",
        "for a layout design; (c) select one period and use the 2-index ",
        "dyadic/missing design."),
        N, names(idx)[1L], dims[1L], names(idx)[2L], dims[2L],
        names(time_v), nT, exp3,
        if (anyDuplicated(tri) > 0L) "; some cells repeat" else ""),
        call. = FALSE)
    exp3
  }

  if (design != "auto") {
    ## forced design: detection is skipped but structural feasibility is
    ## still validated (a forced design that cannot run should fail here,
    ## actionably, not deep inside the engine)
    chosen <- design
    reason <- "forced via design ="
    if (design %in% c("dyadic", "missing")) {
      if (C != 2L) stop(sprintf("design = \"%s\" needs exactly 2 index dimensions.", design), call. = FALSE)
      if (dup2) stop(sprintf(paste0("design = \"%s\" requires one observation per cell, ",
                                    "but some (row, col) cells repeat. Use layout (rep) or panel (time)."),
                             design), call. = FALSE)
      roles <- list(row = names(idx)[1L], col = names(idx)[2L])
      K_default <- min(dims[1:2]) - 1L
      cells_obs <- N; cells_exp <- prod(dims[1:2])
      balance <- if (N == cells_exp) "complete" else sprintf("incomplete (%d of %d cells)", N, cells_exp)
      if (design == "dyadic" && N != cells_exp)
        stop(sprintf(paste0("design = \"dyadic\" requires a complete array but only ",
                            "%d of %d cells are observed. Use design = \"missing\" ",
                            "(Procedure 2, fully observed bicliques) instead."),
                     N, cells_exp), call. = FALSE)
      if (design == "missing") {
        K_default <- NA_integer_          # depends on the biclique blocks
      }
    } else if (design == "panel") {
      if (is.null(time_v) && C != 3L)
        stop("design = \"panel\" needs a time dimension: pass 3 index columns or `time =`.", call. = FALSE)
      if (is.null(time_v)) { time_v <- idx[3L]; idx <- idx[1:2]; dims <- dims[1:2]; dense <- dense[1:2] }
      roles <- list(row = names(idx)[1L], col = names(idx)[2L], time = names(time_v))
      K_default <- min(dims[1:2]) - 1L
      cells_obs <- N; cells_exp <- require_complete_panel()
      balance <- "complete"
    } else if (design == "threeway") {
      if (C != 3L) stop("design = \"threeway\" needs exactly 3 index dimensions.", call. = FALSE)
      cells3 <- .cell_code(cbind(dense[[1L]], dense[[2L]], dense[[3L]]))
      cells_obs <- N; cells_exp <- prod(dims)
      if (anyDuplicated(cells3) > 0L || N != cells_exp)
        stop(sprintf(paste0(
          "design = \"threeway\" needs a complete balanced crossed array: %d ",
          "observations vs %s = %d expected cells%s. Curate a balanced subset ",
          "or use the 2-index dyadic/missing design."),
          N, paste(sprintf("%s=%d", names(dims), dims), collapse = " x "),
          cells_exp,
          if (anyDuplicated(cells3) > 0L) "; some cells repeat" else ""),
          call. = FALSE)
      roles <- list(id1 = names(idx)[1L], id2 = names(idx)[2L], id3 = names(idx)[3L])
      K_default <- min(dims) - 1L
      balance <- "complete"
    } else if (design == "layout") {
      if (C != 2L) stop("design = \"layout\" needs exactly 2 index dimensions (cells).", call. = FALSE)
      finish_layout("forced via design =")
      reason <- "forced via design ="
    }
  } else if (!is.null(rep_v)) {
    ## user declared within-cell replication
    if (C != 2L) stop("With `rep =`, pass exactly 2 index dimensions (the cells).", call. = FALSE)
    finish_layout("`rep =` declares within-cell replication")
  } else if (!is.null(time_v)) {
    ## user declared a time dimension -> panel
    if (C != 2L) stop("With `time =`, pass exactly 2 index dimensions (the cross-section).", call. = FALSE)
    chosen <- "panel"
    roles <- list(row = names(idx)[1L], col = names(idx)[2L], time = names(time_v))
    reason <- "`time =` declares a panel"
    K_default <- min(dims[1:2]) - 1L
    cells_obs <- N; cells_exp <- require_complete_panel()
    balance <- "complete"
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
      cells_obs <- N; cells_exp <- prod(dims[1:2])
      roles <- list(row = names(idx)[1L], col = names(idx)[2L])
      K_default <- min(dims[1:2]) - 1L
      if (N == cells_exp) {
        chosen <- "dyadic"; balance <- "complete"
        reason <- "2 indices, one observation per cell, complete array"
      } else {
        chosen <- "missing"; balance <- sprintf("incomplete (%d of %d cells)", N, cells_exp)
        reason <- sprintf("2 indices, %d of %d cells observed", N, cells_exp)
        K_default <- NA_integer_          # depends on the biclique blocks
        notes <- c(notes, paste0(
          "The permutation-group order under missingness is set by the fully ",
          "observed biclique blocks; see find_bicliques() for the achievable K."))
      }
    }
  } else {
    ## C == 3, no tags: complete balanced crossed array required
    cells3 <- .cell_code(cbind(dense[[1L]], dense[[2L]], dense[[3L]]))
    cells_obs <- N; cells_exp <- prod(dims)
    if (anyDuplicated(cells3) > 0L || N != cells_exp) {
      stop(sprintf(paste0(
        "3 index dimensions but not a complete balanced crossed array ",
        "(%d observations vs %s = %d expected cells%s). mwperm_panel()/",
        "mwperm_threeway() require completeness. Options: (a) curate a ",
        "complete balanced subset; (b) if the third index is replication, ",
        "pass it as `rep =` for a layout design; (c) drop to the 2-index ",
        "dyadic/missing design by aggregating or selecting one level."),
        N, paste(sprintf("%s=%d", names(dims), dims), collapse = " x "),
        cells_exp,
        if (anyDuplicated(cells3) > 0L) "; some cells repeat" else ""),
        call. = FALSE)
    }
    ## which index is time? tagged > name > class/values > ambiguous
    name_hit <- which(vapply(names(idx), .timelike_name, logical(1)))
    val_hit <- which(vapply(seq_len(3L), function(k)
      .timelike_values(idx[[k]], dims[k], max(dims[-k])), logical(1)))
    t_k <- NULL; t_why <- NULL
    if (length(name_hit) == 1L) {
      t_k <- name_hit; t_why <- sprintf("'%s' identified as time by name", names(idx)[t_k])
    } else if (length(val_hit) == 1L) {
      t_k <- val_hit
      t_why <- sprintf("'%s' identified as time by its values (temporal/regularly spaced)",
                       names(idx)[t_k])
    }
    if (is.null(t_k)) {
      ## ASSUMPTION FORK: ambiguous -> panel default on the third index
      t_k <- 3L
      t_why <- sprintf("ambiguous third dimension; '%s' held fixed by default", names(idx)[3L])
      warns <- c(warns, paste0(
        "No index is clearly the time dimension, so the PANEL design was ",
        "chosen with '", names(idx)[3L], "' held fixed (valid even if all ",
        "three dimensions are exchangeable, only less powerful). If all ",
        "three really are exchangeable, force design = \"threeway\"; if a ",
        "different index is time, tag it via time = ."))
    }
    chosen <- "panel"
    time_v <- idx[t_k]; idx <- idx[-t_k]
    dims_cs <- dims[-t_k]
    roles <- list(row = names(idx)[1L], col = names(idx)[2L], time = names(time_v))
    reason <- t_why
    K_default <- min(dims_cs) - 1L
    balance <- "complete"
    dims <- c(dims_cs, dims[t_k])
    names(dims) <- c(names(idx), names(time_v))
  }

  K_default <- if (is.na(K_default)) NA_integer_ else min(K_default, 199L)
  res_ok <- if (is.na(K_default)) NA else (K_default + 1L) >= 20L

  ## the exact downstream call
  fn <- c(dyadic = "mwperm_dyadic", panel = "mwperm_panel",
          threeway = "mwperm_threeway", layout = "mwperm_layout",
          missing = "mwperm_missing")[[chosen]]
  args <- switch(chosen,
    dyadic = sprintf("row = %s, col = %s", roles$row, roles$col),
    missing = sprintf("row = %s, col = %s, min_block = ...", roles$row, roles$col),
    panel = sprintf("row = %s, col = %s, time = %s, time_fe = TRUE",
                    roles$row, roles$col, roles$time),
    threeway = sprintf("id1 = %s, id2 = %s, id3 = %s",
                       roles$id1, roles$id2, roles$id3),
    layout = sprintf("row = %s, col = %s%s", roles$row, roles$col,
                     if (identical(roles$rep, "(within-cell order)")) ""
                     else sprintf(", rep = %s", roles$rep)))
  call_str <- sprintf("%s(y, d, x, %s)", fn, args)

  structure(list(
    design = chosen, roles = roles, dims = dims, n_obs = N,
    cells = c(observed = cells_obs, expected = cells_exp),
    balance = balance, K_default = K_default, resolution_ok = res_ok,
    call_str = call_str, reason = reason,
    warnings = warns, notes = notes,
    time = if (!is.null(time_v)) time_v[[1L]] else NULL,
    rep = if (!is.null(rep_v)) rep_v[[1L]] else NULL,
    index = idx
  ), class = "mwperm_design")
}

#' Print a design diagnosis
#'
#' @param x An object of class \code{"mwperm_design"} from
#'   \code{\link{mwperm_check}}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.mwperm_design <- function(x, ...) {
  cat("\nmwperm design diagnosis\n")
  cat(strrep("-", 30), "\n", sep = "")
  cat("Detected design :", x$design, sprintf("(%s)\n", x$reason))
  role_txt <- paste(sprintf("%s = %s", names(x$roles), unlist(x$roles)),
                    collapse = ", ")
  cat("Roles           :", role_txt, "\n")
  cat("Dimensions      :",
      paste(sprintf("%s (%d)", names(x$dims), x$dims), collapse = " x "),
      sprintf("| %d observations\n", x$n_obs))
  cat("Balance         :", x$balance, "\n")
  if (is.na(x$K_default)) {
    cat("Resolution      : depends on the biclique blocks (see find_bicliques)\n")
  } else {
    cat(sprintf("Resolution      : default K = %d, smallest attainable p = 1/%d %s\n",
                x$K_default, x$K_default + 1L,
                if (isTRUE(x$resolution_ok)) "(95% confidence set attainable)"
                else "(TOO COARSE for a 95% confidence set; needs >= 20 clusters)"))
  }
  cat("Would run       :", x$call_str, "\n")
  for (w in x$warnings)
    cat(strwrap(w, initial = "  ! ", prefix = "    ",
                width = 0.9 * getOption("width", 80)), sep = "\n")
  for (nt in x$notes)
    cat(strwrap(nt, initial = "  - ", prefix = "    ",
                width = 0.9 * getOption("width", 80)), sep = "\n")
  cat("\n")
  invisible(x)
}

## ---- mwperm: one-call dispatch ----------------------------------------------

#' One-call invariant permutation test with automatic design detection
#'
#' Detects the clustering design of the data via \code{\link{mwperm_check}}
#' and dispatches to the matching test -- \code{\link{mwperm_dyadic}},
#' \code{\link{mwperm_panel}}, \code{\link{mwperm_threeway}},
#' \code{\link{mwperm_layout}} or \code{\link{mwperm_missing}} -- forwarding
#' all arguments unchanged. A thin convenience layer: the returned object is
#' exactly what the underlying function returns (plus a record of what was
#' detected), and calling the specific function directly with the same seed
#' gives identical results.
#'
#' See \code{\link{mwperm_check}} for the detection rules, in particular the
#' two assumption-dependent forks (panel-vs-threeway and layout-vs-panel)
#' that are announced rather than silently resolved.
#'
#' @param y,d,x Outcome, covariate(s) of interest, and optional nuisance
#'   covariates, as in \code{\link{mwperm_dyadic}}. With \code{data} given,
#'   each may also be a character (vector of) column name(s) resolved
#'   against it.
#' @param index The clustering dimensions (2 or 3): a data frame, named list
#'   of vectors, or character vector of column names in \code{data}.
#' @param data Optional data frame; column names in \code{y}, \code{d},
#'   \code{x}, \code{index}, \code{time}, \code{rep} are resolved against it.
#' @param time,rep Optional explicit role tags (vector or column name); see
#'   \code{\link{mwperm_check}}.
#' @param design Force a design instead of auto-detecting.
#' @param time_fe Passed to \code{\link{mwperm_panel}} (panel only).
#' @param L0 Passed to \code{\link{mwperm_layout}} (layout only).
#' @param min_block,block_method Passed to \code{\link{mwperm_missing}}
#'   (missing only).
#' @param verbose If \code{TRUE} (default) print one line stating the
#'   detected design and the dispatched call.
#' @inheritParams mwperm_dyadic
#'
#' @return The \code{"mwperm"} object of the dispatched test, with an extra
#'   \code{auto} field recording the detection (design, roles, reason); the
#'   detection notices are prepended to the object's \code{note} field and
#'   shown by \code{\link{print.mwperm}}.
#'
#' @examples
#' data(trade_dyadic)
#' fit <- mwperm(y = "log_trade", d = "log_dist",
#'               x = c("log_gdp_i", "log_gdp_j"),
#'               index = c("importer", "exporter"),
#'               data = trade_dyadic, seed = 1)
#' fit
#' @export
mwperm <- function(y, d, x = NULL, index, data = NULL, time = NULL, rep = NULL,
                   design = c("auto", "dyadic", "threeway", "panel", "layout",
                              "missing"),
                   K = NULL, alpha = 0.05, beta_null = 0, conf_int = TRUE,
                   n_reps = 1L, seed = NULL, grid = NULL, n_cores = 1L,
                   time_fe = TRUE, L0 = NULL, min_block = 3L,
                   block_method = c("greedy", "exact"), verbose = TRUE) {
  design <- match.arg(design)
  cl <- match.call()
  ## capture the caller's expression for d BEFORE evaluation: the front ends
  ## label coefficients by deparse(substitute(d)), which through do.call would
  ## deparse the VALUES; we forward d as a named-column matrix instead.
  d_expr <- paste(deparse(substitute(d)), collapse = "")

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
      stop(sprintf(paste0("`%s` is a character vector; to use column names, pass ",
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
                      time = time, rep = rep, design = design)

  ## design-specific arguments must not be silently accepted for the wrong
  ## design (warn and ignore, mirroring the printed diagnosis)
  supplied <- names(cl)
  check_arg <- function(arg, ok_design) {
    if (arg %in% supplied && chk$design != ok_design)
      warning(sprintf("`%s` applies to the %s design only; it was ignored for '%s'.",
                      arg, ok_design, chk$design), call. = FALSE)
  }
  check_arg("time_fe", "panel")
  check_arg("L0", "layout")
  check_arg("min_block", "missing")
  check_arg("block_method", "missing")

  if (isTRUE(verbose)) {
    message(sprintf("Detected design: %s (%s) -> running %s",
                    chk$design, chk$reason, chk$call_str))
  }

  common <- list(y = y, d = d, x = x, K = K, alpha = alpha,
                 beta_null = beta_null, conf_int = conf_int, n_reps = n_reps,
                 seed = seed, grid = grid, n_cores = n_cores)
  ix <- chk$index
  res <- switch(chk$design,
    dyadic = do.call(mwperm_dyadic,
                     c(common, list(row = ix[[1L]], col = ix[[2L]]))),
    missing = do.call(mwperm_missing,
                      c(common, list(row = ix[[1L]], col = ix[[2L]],
                                     min_block = min_block,
                                     block_method = block_method))),
    panel = do.call(mwperm_panel,
                    c(common, list(row = ix[[1L]], col = ix[[2L]],
                                   time = chk$time, time_fe = time_fe))),
    threeway = do.call(mwperm_threeway,
                       c(common, list(id1 = ix[[1L]], id2 = ix[[2L]],
                                      id3 = ix[[3L]]))),
    layout = do.call(mwperm_layout,
                     c(common, list(row = ix[[1L]], col = ix[[2L]],
                                    rep = chk$rep, L0 = L0))))

  ## make the automatic choice transparent on the returned object
  res$auto <- list(design = chk$design, reason = chk$reason, roles = chk$roles)
  res$note <- c(chk$warnings, chk$notes, res$note)
  res$call <- cl
  res
}
