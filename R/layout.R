#' Invariant permutation test for two-way layouts (within-cell replication)
#'
#' Finite-sample valid test of \eqn{H_0: \beta = b} for a two-way layout
#' \deqn{y_{ijl} = x_{ijl}^\top \gamma + d_{ijl}^\top \beta + \varepsilon_{ijl}}
#' in which each \eqn{(i, j)} cell contains \eqn{\ell_{ij}} replicate
#' observations and the cell sizes may be \emph{unequal}. Permutations are
#' applied only within cells, over the replication index \eqn{l}. This is valid
#' under the relaxed condition that errors are exchangeable with respect to
#' \eqn{l} within each cell (condition InvA restricted to \eqn{l}), e.g.
#' \eqn{\varepsilon_{ijl} = \eta_{ij} + \zeta_l + u_{ijl}} with arbitrary
#' \eqn{\eta_{ij}} but i.i.d. \eqn{\zeta_l}. It is appropriate when \eqn{l}
#' indexes independent replications (two-way layouts in randomised
#' experiments). See Guo, Toulis and Wang (2026), Section 6.3.
#'
#' Note: the covariate of interest must vary within cells. If \eqn{d_{ijl}} is
#' constant within every \eqn{(i, j)} cell (a cell-level covariate), within-cell
#' permutation yields a trivial test with no power; a warning is issued.
#'
#' @inheritParams mwperm_dyadic
#' @param row,col Cell identifiers along the two layout dimensions.
#' @param rep Optional replication identifier within each cell; if \code{NULL},
#'   the order of appearance within a cell is used.
#' @param L0 Optional integer capacity threshold for balancing an unbalanced
#'   layout (Section 6.3 of Guo, Toulis and Wang, 2026). When supplied, cells
#'   with fewer than \code{L0} replicates are dropped and each remaining cell is
#'   uniformly downsampled to exactly \code{L0} replicates, giving a balanced
#'   array before testing. The downsampling is reproducible through \code{seed}.
#'   When \code{NULL} (the default) all replicates are used with their (possibly
#'   unequal) cell sizes.
#' @param K Number of non-identity permutations; defaults to
#'   \code{min(cell size) - 1} capped at 199 (with \code{L0}, to \code{L0 - 1}).
#'   Must satisfy \code{K + 1 <= min(cell size)}.
#'
#' @return An object of class \code{"mwperm"}: \code{estimate}/\code{se_naive}
#'   are the OLS estimate and naive SE, \code{conf_int} (or
#'   \code{conf_region}/\code{conf_box} for several coefficients) the IPT
#'   inverted-test confidence set, and \code{pvalue} the IPT permutation
#'   p-value; see \code{\link{mwperm_dyadic}} for the field provenance in full.
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation
#'   inference under multi-way clustering and missing data, Section 6.3.
#'   arXiv:2601.08610.
#' @seealso \code{\link{mwperm_dyadic}}, \code{\link{mwperm_missing}}.
#' @export
mwperm_layout <- function(y, d, x = NULL, row, col, rep = NULL, L0 = NULL,
                          K = NULL, alpha = 0.05, beta_null = 0, conf_int = TRUE,
                          n_reps = 1L, seed = NULL, grid = NULL, n_cores = 1L) {
  cl <- match.call()
  y <- .check_y(y); N <- length(y)
  D <- as.matrix(d); d_names <- .coef_names(D, deparse(substitute(d)))
  X <- .make_X(x, N)
  .check_lengths(N, list(row = row, col = col,
                         rep = if (is.null(rep)) seq_len(N) else rep))

  cell <- .dense_id(interaction(.dense_id(row, "row"), .dense_id(col, "col"),
                                drop = TRUE))  # dense (row,col) cell id
  ncell <- max(cell)                   # number of occupied cells
  ell <- tabulate(cell, nbins = ncell) # replicate count per cell (the cell sizes)

  ## optional balancing: keep dense cells, downsample each to exactly L0 (5.3)
  balance_note <- character(0)
  if (!is.null(L0)) {
    L0 <- as.integer(L0)
    if (length(L0) != 1L || is.na(L0) || L0 < 2L)
      stop("`L0` must be a single integer >= 2.", call. = FALSE)
    if (!any(ell >= L0))
      stop(sprintf("No cell has at least L0 = %d replicates (largest cell has %d).",
                   L0, max(ell)), call. = FALSE)
    keep <- .downsample_to_L0(cell, ell, L0, seed = seed)
    balance_note <- sprintf(paste0(
      "Balanced to L0 = %d replicates/cell: kept %d of %d cells and %d of %d ",
      "observations (cells with < L0 replicates dropped; dense cells uniformly ",
      "downsampled)."), L0, sum(ell >= L0), ncell, length(keep), N)
    y <- y[keep]; D <- D[keep, , drop = FALSE]; X <- X[keep, , drop = FALSE]
    if (!is.null(rep)) rep <- rep[keep]
    N <- length(y)
    cell <- .dense_id(cell[keep]); ncell <- max(cell)
    ell <- tabulate(cell, nbins = ncell)
  }

  ## Validate here (the engine re-checks) because the within-cell-variation
  ## diagnostic below reads `d` and would fail opaquely on NA/non-numeric input.
  .check_finite(list(y = y, d = D, x = X))

  ## Dense 1-based replication index l within each cell. The permutation acts on
  ## this index, so it must run 1..ell[c] inside every cell; ties in the ordering
  ## key are broken deterministically.
  widx <- integer(N)
  ord_key <- if (is.null(rep)) seq_len(N) else as.numeric(factor(rep))  # within-cell ordering
  for (c in seq_len(ncell)) {
    idx <- which(cell == c)
    widx[idx] <- rank(ord_key[idx], ties.method = "first")
  }

  ## Warn if d does not vary within any cell: within-cell permutation then leaves
  ## the residual statistic unchanged, so the test would have no power.
  within_var <- tapply(seq_len(N), cell, function(ii)
    any(apply(D[ii, , drop = FALSE], 2L, function(v) diff(range(v)) > 0)))
  if (!any(unlist(within_var)))
    warning("`d` is constant within every cell; the within-cell test has no power. ",
            "Did you mean mwperm_dyadic() or mwperm_panel()?", call. = FALSE)

  min_cell <- min(ell)
  K <- .default_K(K, min_cell)

  perm_builder <- function(rep_seed) {
    cell_groups <- vector("list", ncell)
    for (c in seq_len(ncell))
      cell_groups[[c]] <- build_perm_set(ell[c], K, seed = .sub_seed(rep_seed, c))
    .build_obs_perms_layout(cell, widx, cell_groups)
  }

  res <- .ipt_engine(y, D, X, perm_builder, K = K, n_reps = n_reps, seed = seed,
                     alpha = alpha, conf_int = conf_int, beta_null = beta_null,
                     grid = grid,
                     type = if (is.null(L0)) "layout" else "layout (balanced)",
                     d_names = d_names,
                     n_clusters = c(ncell = ncell, min_cell = min_cell), call = cl,
                     n_cores = n_cores)
  if (length(balance_note)) res$note <- c(balance_note, res$note)
  res
}

#' Choose row indices balancing every dense cell to exactly L0 replicates.
#'
#' Keeps cells with at least \code{L0} observations and, within each, draws a
#' uniform random subset of size \code{L0}. Reproducible via \code{seed} with
#' the usual RNG-state hygiene (the global stream is left untouched).
#' @keywords internal
#' @noRd
.downsample_to_L0 <- function(cell, ell, L0, seed = NULL) {
  if (!is.null(seed)) {
    old <- .save_seed()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }
  dense <- which(ell >= L0)            # cells big enough to keep
  per_cell <- split(seq_along(cell), cell)   # observation indices grouped by cell
  keep <- vector("list", length(dense))
  for (q in seq_along(dense)) {
    idx <- per_cell[[dense[q]]]        # all rows of this dense cell
    keep[[q]] <- idx[sample.int(length(idx), L0)]   # uniform subsample of size L0
  }
  sort(unlist(keep, use.names = FALSE))
}

#' Within-cell observation permutations for two-way layouts.
#' @keywords internal
#' @noRd
.build_obs_perms_layout <- function(cell, widx, cell_groups) {
  N <- length(cell)
  ncell <- max(cell)
  Kp1 <- length(cell_groups[[1L]])     # group order (all cells share it)
  cell_idx <- split(seq_len(N), cell)            # global observation positions per cell
  ## pos_in_cell[[c]][l] = the global observation index whose within-cell index is l,
  ## i.e. a lookup from replication index back to row, per cell.
  pos_in_cell <- lapply(seq_len(ncell), function(c) {
    idx <- cell_idx[[as.character(c)]]
    idx[order(widx[idx])]
  })
  obs_perms <- vector("list", Kp1)
  for (k in seq_len(Kp1)) {
    g <- integer(N)                    # observation gather-vector for element k
    for (c in seq_len(ncell)) {
      idx <- cell_idx[[as.character(c)]]
      img_within <- cell_groups[[c]][[k]][widx[idx]]   # permuted within-cell indices
      g[idx] <- pos_in_cell[[c]][img_within]           # back to global observation indices
    }
    obs_perms[[k]] <- g
  }
  obs_perms
}
