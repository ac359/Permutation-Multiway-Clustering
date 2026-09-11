#' Invariant permutation test for incomplete panels
#'
#' Finite-sample valid test of H0: beta = b for a panel
#' \preformatted{  y_ijt = x_ijt' gamma + d_ijt' beta + eps_ijt}
#'
#' in which the (i, j) array is **incomplete** -- some pairs are never
#' observed, or are observed in only some periods -- so [mwperm_panel()], which
#' requires a complete balanced array, cannot run. This combines condition
#' InvB of Guo, Toulis and Wang (2026), Section 6.2 with the missing-data
#' machinery of Section 5, which is the extension their Section 9 leaves open.
#'
#' The construction is:
#' 1. form the mask `M_ij = 1` if pair (i, j) is observed in EVERY period, and
#'   0 otherwise;
#' 2. run the biclique search ([find_bicliques()], Algorithm 2) on that mask to
#'   obtain disjoint fully observed blocks `I_q` by `J_q`;
#' 3. discard every observation outside the selected blocks;
#' 4. apply Procedure 2 to what remains: a row group on `I_q` and a column
#'   group on `J_q` per block, of common order K + 1, applied identically in
#'   every period, so cell (i, j) in period t maps to cell (pi(i), sigma(j)) in
#'   period t.
#'
#' Holding the period fixed is the same device [mwperm_panel()] uses, and it is
#' what keeps the test valid when the errors are autocorrelated over time.
#'
#' @section Assumptions: Two conditions, and neither is within-period
#'   independence. The first is condition InvB applied blockwise: for every
#'   period t, and with the *same* (pi, sigma) used in every period,
#'
#' (eps_ijt for i in I_q, j in J_q) has the same distribution, given X and D,
#' as (eps_[pi(i)][sigma(j)]t for i in I_q, j in J_q)
#'
#'   The second is Assumption 4 on the mask: which pairs are completely
#'   observed is independent of the errors given the covariates. The mask may
#'   depend on the covariates in any way -- an unbalanced panel where richer
#'   country pairs are observed in more years is fine -- but not on the
#'   outcome. Under those two conditions the p-value is exact in finite
#'   samples, by the argument of Theorem 4 applied within each block.
#'
#'   Arbitrary common time effects are permitted, exactly as in
#'   [mwperm_panel()]: nothing is assumed across periods, so a trend, a break,
#'   or serial correlation of any form costs no validity.
#'
#' @section What is discarded: A pair observed in some periods but not all
#'   fails the mask and is dropped whole, and pairs that clear the mask are
#'   still dropped if the biclique search cannot place them in a block. That
#'   loss is the price of exactness under missingness, and it is reported in
#'   `note`, `cells_used` and `cells_total`. Dropping *periods* instead --
#'   keeping only the years in which every pair is observed -- is often the
#'   better trade and is done by subsetting the data before calling this
#'   function; there is no free lunch, and which side to cut depends on the
#'   panel.
#'
#' @param y Numeric outcome vector.
#' @param d Numeric vector or matrix of covariate(s) of interest. With a single
#'   covariate a confidence interval is returned; with several, a joint
#'   acceptance region.
#' @param x Optional numeric matrix or data frame of nuisance covariates. An
#'   intercept is added automatically.
#' @param row,col Cluster identifiers of the two cross-sectional dimensions,
#'   for example importer and exporter. Any type; recoded to dense ids.
#' @param time Period identifier. Held fixed by the permutation -- the errors
#'   are never assumed exchangeable across it.
#' @param K Number of non-identity permutations; defaults to `min(smallest
#'   block side) - 1` capped at 199. Must satisfy `K + 1 <=` the smallest side
#'   of the smallest selected block.
#' @param min_block Minimum block side for the biclique search, a single
#'   integer or a length-2 integer giving the row and column floors
#'   separately. Floored at 2.
#' @param block_method `"greedy"` (default) or `"exact"`; see
#'   [find_bicliques()]. A sub-maximal block costs power, never validity.
#' @param alpha Significance level used for the confidence set and the printed
#'   decision.
#' @param beta_null Null value(s); length 1 or `ncol(d)`.
#' @param conf_int Logical; compute the confidence set by test inversion.
#' @param n_reps Number of independent permutation groups; the reported
#'   p-value is aggregated over them.
#' @param seed Integer seed. Required for reproducibility, and for the
#'   rep-parallel path.
#' @param grid Optional explicit grid of null values for the inversion.
#' @param time_fe Logical; if `TRUE` (the default) period dummies are added to
#'   the nuisance design. They are invariant to the within-period permutation,
#'   so they cost no validity, and they remove the common time effect from the
#'   point estimate.
#' @param aggregate Cross-repetition rule, `"median"` (default) or
#'   `"median2"`.
#' @param n_cores Integer; 1 (default) runs serially and is bit-identical to
#'   any larger value.
#'
#' @return An object of class `"mwperm"`. Beyond the usual fields it carries
#'   `n_blocks`, `cells_used` and `cells_total`, the block accounting.
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data. Condition InvB is Section
#'   6.2, Procedure 2 and Assumption 4 are Section 5. arXiv:2601.08610.
#'
#' @seealso [mwperm_panel()] for a complete balanced panel, [mwperm_missing()]
#'   for an incomplete array with no time dimension, and [find_bicliques()]
#'   for the block search.
#'
#' @examples
#' data(trade_panel)
#' ## thin 20 country pairs so they are observed in the first year only:
#' ## an unbalanced panel, which mwperm_panel() cannot test
#' set.seed(1)
#' tp <- trade_panel
#' pair <- paste(tp$importer, tp$exporter)
#' thin <- sample(unique(pair), 20)
#' tp <- tp[!(pair %in% thin & tp$year > min(tp$year)), ]
#' fit <- with(tp, mwperm_panel_missing(y = log_trade, d = fta,
#'                                      x = cbind(log_gdp_i, log_gdp_j),
#'                                      row = importer, col = exporter,
#'                                      time = year, min_block = 5, seed = 1))
#' fit
#' ## Note what the fit reports: 20 thinned pairs cost far more than 20 cells,
#' ## because a block must be complete in BOTH margins. With 22 clusters per
#' ## side an incomplete panel rarely supports a 95% interval -- K + 1 >= 20
#' ## needs a fully observed 20 x 20 block. Dropping the sparsest periods
#' ## before calling is usually the better trade.
#' @export
mwperm_panel_missing <- function(y, d, x = NULL, row, col, time, K = NULL,
                                 min_block = 3L,
                                 block_method = c("greedy", "exact"),
                                 alpha = 0.05, beta_null = 0, conf_int = TRUE,
                                 n_reps = 10L, seed = NULL, grid = NULL,
                                 time_fe = TRUE,
                                 aggregate = c("median", "median2"),
                                 n_cores = 1L) {
  cl <- match.call()
  aggregate <- match.arg(aggregate)
  block_method <- match.arg(block_method)
  y <- .check_y(y)
  N <- length(y)
  D <- as.matrix(d)
  d_names <- .coef_names(D, deparse(substitute(d)))
  X <- .make_X(x, N)
  .check_lengths(N, list(row = row, col = col, time = time))

  ri <- .dense_id(row, "row")
  ci <- .dense_id(col, "col")
  ti <- .dense_id(time, "time")
  n_row <- max(ri)
  n_col <- max(ci)
  n_t <- max(ti)

  ## One observation per (row, col, period). Without this the "same (pi, sigma)
  ## in every period" map is not well defined: two observations in one cell and
  ## period would both claim the same slot, and the gather vector would resolve
  ## the collision silently to whichever came first.
  if (anyDuplicated(cbind(ri, ci, ti)))
    stop(paste0("Each (row, col, time) cell must appear at most once. For ",
                "several observations per cell use mwperm_layout() (if the ",
                "repeats are exchangeable replicates) or mwperm_irregular() ",
                "(if they are not)."), call. = FALSE)

  ## --- (i) the mask: pairs observed in EVERY period -------------------------
  ## Given uniqueness above, a cell has all n_t periods exactly when it holds
  ## n_t observations, so the count is the set condition.
  cell <- .dense_id(interaction(ri, ci, drop = TRUE))
  ncell <- max(cell)
  ell <- tabulate(cell, nbins = ncell)
  first <- match(seq_len(ncell), cell)
  mask <- ell == n_t
  if (!any(mask))
    stop(sprintf(paste0("No (row, col) pair is observed in all %d periods, ",
                        "so the mask is empty (the fullest pair has %d of ",
                        "%d). Drop the sparsest periods and call again, or ",
                        "use mwperm_missing() on a single period."),
                 n_t, max(ell), n_t), call. = FALSE)

  ## --- (ii) biclique search on the mask (Algorithm 2) -----------------------
  blocks <- find_bicliques(ri[first[mask]], ci[first[mask]],
                           min_block = min_block, method = block_method)
  if (length(blocks) == 0L)
    stop(sprintf(paste0("No fully observed block with both sides >= ",
                        "min_block = %s was found among the %d pairs ",
                        "observed in all %d periods. Lower `min_block`, or ",
                        "drop the sparsest periods so that more pairs clear ",
                        "the mask."),
                 paste(min_block, collapse = "x"), sum(mask), n_t),
         call. = FALSE)

  ## --- restrict to the cells inside the selected blocks ---------------------
  in_block <- integer(ncell)              # block index q per cell; 0 = none
  for (q in seq_along(blocks)) {
    b <- blocks[[q]]
    hit <- mask & (ri[first] %in% b$rows) & (ci[first] %in% b$cols)
    in_block[hit] <- q
  }
  idx <- which(in_block[cell] > 0L)       # retained observations
  if (!length(idx))
    stop("Internal error: the selected blocks contain no observations.",
         call. = FALSE)   # nocov

  ## Period dummies are built on the full data and then subset: every period
  ## survives by construction (a retained cell has all of them), so the columns
  ## are the same either way, and this keeps the reference level stable.
  if (isTRUE(time_fe) && n_t > 1L) {
    TD <- stats::model.matrix(~ factor(ti))[, -1L, drop = FALSE]
    colnames(TD) <- paste0("time", sort(unique(ti))[-1L])
    X <- cbind(X, TD)
  }

  yk <- y[idx]
  Dk <- D[idx, , drop = FALSE]
  Xk <- X[idx, , drop = FALSE]
  ri_k <- ri[idx]
  ci_k <- ci[idx]
  ti_k <- ti[idx]                          # the slot, held fixed
  Nk <- length(idx)
  n_cells_used <- sum(in_block > 0L)

  ## Block index and block-local row/column position of every retained
  ## observation, for the Procedure 2 permutation builder.
  blk_k <- in_block[cell[idx]]
  lrow_k <- integer(Nk)
  lcol_k <- integer(Nk)
  for (q in seq_along(blocks)) {
    s <- which(blk_k == q)
    if (!length(s)) next
    lrow_k[s] <- match(ri_k[s], blocks[[q]]$rows)
    lcol_k[s] <- match(ci_k[s], blocks[[q]]$cols)
  }

  ## --- group order: the smallest permuted block side over the blocks --------
  min_side <- min(vapply(blocks,
                         function(b) min(length(b$rows), length(b$cols)),
                         integer(1)))
  K_was_null <- is.null(K)
  K <- .default_K(K, min_side)

  sizes <- vapply(blocks, function(b) c(length(b$rows), length(b$cols)),
                  integer(2))
  note <- c(
    sprintf(paste0("Incomplete panel: %d of %d observed (row, col) pairs are ",
                   "present in all %d periods, and %d of those lie in the %d ",
                   "fully observed block%s the biclique search extracted. ",
                   "Kept %d of %d observations (%.1f%%). The discarded pairs ",
                   "are what buys exact validity under missingness."),
            sum(mask), ncell, n_t, n_cells_used, length(blocks),
            if (length(blocks) == 1L) "" else "s", Nk, N, 100 * Nk / N),
    sprintf("Block sizes (rows x cols): %s.",
            paste(sprintf("%dx%d", sizes[1, ], sizes[2, ]), collapse = ", ")),
    if (K_was_null && is.numeric(alpha) && length(alpha) == 1L &&
        is.finite(alpha) && alpha > 0 && alpha < 1 && 1 / (K + 1) > alpha)
      sprintf(paste0("Resolution here is set by the smallest selected ",
                     "block: its permuted side is %d, so K = %d. Raise ",
                     "`min_block` so that small blocks cannot set K, or drop ",
                     "the sparsest periods so that more pairs clear the ",
                     "mask."),
              min_side, K)
    else character(0))

  ## --- (iv) Procedure 2 with the period held fixed --------------------------
  perm_builder <- function(rep_seed)
    .build_obs_perms_blocks(rep_seed, K, blocks,
                            ri = ri_k, ci = ci_k, blk = blk_k,
                            lrow = lrow_k, lcol = lcol_k, permute = "both",
                            slot = ti_k)

  res <- .ipt_engine(yk, Dk, Xk, perm_builder, K = K, n_reps = n_reps,
                     seed = seed, alpha = alpha, conf_int = conf_int,
                     beta_null = beta_null, grid = grid,
                     type = "panel (bicliques)", d_names = d_names,
                     n_clusters = c(row = n_row, col = n_col, time = n_t),
                     call = cl, n_cores = n_cores, ci_agg = aggregate)
  res$note <- c(note, res$note)
  res$n_blocks <- length(blocks)
  res$cells_used <- n_cells_used
  res$cells_total <- ncell
  res
}
