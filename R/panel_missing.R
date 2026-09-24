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
#'   0 otherwise (with `L0`, in every period of a common set S of L0 periods;
#'   see *Keeping L0 periods*);
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
#' @section Keeping L0 periods (`L0`): A pair that misses even one period
#'   fails the every-period mask. When many pairs miss a few periods it is
#'   often the better trade to keep fewer periods and more pairs: with an
#'   integer `L0` the mask becomes "pair observed in every period of S",
#'   where S is the set of L0 periods jointly observed by the most pairs,
#'   chosen greedily from the observation pattern alone (the period observed
#'   by the most pairs first, then the period that keeps the most pairs
#'   observing the whole set; ties to the earlier period), and observations
#'   outside S are dropped. S is a function of the mask and never of the
#'   outcome, which is what Assumption 4 needs, so the cut costs no
#'   validity; the period dummies are built on S. Period s is the *same*
#'   period in every cell, and the permutation holds it fixed -- exactly
#'   condition InvB. `L0` equal to the number of periods reproduces the
#'   default; the fit reports S in `periods_used` and in `note`.
#'
#'   This is also the incomplete-panel reading of the irregular design of
#'   Section 6.4 -- one observation per cell and period -- and it is the
#'   route for any repeats that are *periods*: [mwperm_irregular()]'s random
#'   per-cell trim keeps different periods in different cells, is exact
#'   only for exchangeable replicates, and can over-reject under a period
#'   effect even with a cell-constant `d` (0.39 at nominal 0.05 with
#'   staggered observation windows, against 0.028 here; 1.000 against 0.040
#'   when `d` also varies over the periods). Until 0.4.1
#'   this cut was `mwperm_irregular(trim = "levels")`;
#'   `mwperm_panel_missing(time = <rep>, L0 = <L0>, time_fe = FALSE)`
#'   reproduces it exactly.
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
#'   keeping the L0 years observed by the most pairs -- is often the better
#'   trade and is what `L0 =` does; there is no free lunch, and which side to
#'   cut depends on the panel.
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
#' @param L0 `NULL` (the default) keeps the pairs observed in *every* period.
#'   An integer between 2 and the number of periods keeps instead the L0
#'   periods jointly observed by the most pairs, the same ones in every cell,
#'   and the pairs observed in all of them; see *Keeping L0 periods*.
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
#'   `n_blocks`, `cells_used` and `cells_total`, the block accounting, and
#'   `periods_used`, the labels of the L0 retained periods (`NULL` when `L0`
#'   is `NULL`).
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data. Condition InvB is Section
#'   6.2, Procedure 2 and Assumption 4 are Section 5. arXiv:2601.08610.
#'
#' @seealso [mwperm_panel()] for a complete balanced panel, [mwperm_missing()]
#'   for an incomplete array with no time dimension, [mwperm_irregular()] for
#'   several exchangeable replicates per cell, and [find_bicliques()] for the
#'   block search.
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
#' ## is usually the better trade, and `L0 =` does it: when the thinned pairs
#' ## miss only the LAST year, keeping the 5 best-covered years keeps them all
#' tp5 <- trade_panel[!(pair %in% thin &
#'                        trade_panel$year == max(trade_panel$year)), ]
#' fit5 <- with(tp5, mwperm_panel_missing(y = log_trade, d = fta,
#'                                        x = cbind(log_gdp_i, log_gdp_j),
#'                                        row = importer, col = exporter,
#'                                        time = year, L0 = 5, min_block = 5,
#'                                        seed = 1))
#' fit5$periods_used
#' @export
mwperm_panel_missing <- function(y, d, x = NULL, row, col, time, L0 = NULL,
                                 K = NULL, min_block = 3L,
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
  if (!is.null(L0)) {
    L0 <- suppressWarnings(as.integer(L0))
    if (length(L0) != 1L || is.na(L0) || L0 < 2L)
      stop("`L0` must be a single integer >= 2.", call. = FALSE)
  }

  ## --- (i)-(iii) the mask, the blocks and the retained observations --------
  pd <- .panel_missing_design(row, col, time, L0 = L0, min_block = min_block,
                              block_method = block_method)
  idx <- pd$idx
  blocks <- pd$blocks
  n_t <- pd$n_t
  S <- pd$S                                 # retained period ids (dense)

  ## Period dummies on the retained periods S (every retained cell observes
  ## all of them, so the columns are the same whichever cell defines them and
  ## the reference level, the lowest retained period, is stable). With
  ## L0 = NULL, S is every period and this is the historical design.
  Xk <- X[idx, , drop = FALSE]
  if (isTRUE(time_fe) && length(S) > 1L) {
    TD <- stats::model.matrix(~ factor(pd$ti_k, levels = S))[, -1L,
                                                            drop = FALSE]
    colnames(TD) <- paste0("time", S[-1L])
    Xk <- cbind(Xk, TD)
  }
  yk <- y[idx]
  Dk <- D[idx, , drop = FALSE]
  Nk <- length(idx)
  n_cells_used <- pd$n_cells_used

  ## --- group order: the smallest permuted block side over the blocks --------
  K_was_null <- is.null(K)
  K <- .default_K(K, pd$min_side)

  sizes <- vapply(blocks, function(b) c(length(b$rows), length(b$cols)),
                  integer(2))
  note <- c(
    if (is.null(L0))
      sprintf(paste0("Incomplete panel: %d of %d observed (row, col) pairs ",
                     "are present in all %d periods, and %d of those lie in ",
                     "the %d fully observed block%s the biclique search ",
                     "extracted. Kept %d of %d observations (%.1f%%). The ",
                     "discarded pairs are what buys exact validity under ",
                     "missingness."),
              pd$n_mask, pd$ncell, n_t, n_cells_used, length(blocks),
              if (length(blocks) == 1L) "" else "s", Nk, N, 100 * Nk / N)
    else {
      lev_txt <- if (L0 <= 10L) paste(pd$S_labels, collapse = ", ")
                 else sprintf("%s, ..., %s",
                              paste(pd$S_labels[1:3], collapse = ", "),
                              pd$S_labels[L0])
      sprintf(paste0("Incomplete panel, L0 = %d: the retained periods are ",
                     "{%s}, the %d of the %d periods jointly observed by the ",
                     "most (row, col) pairs; %d of %d observed pairs are ",
                     "present in all of them, %d of those lie in the %d fully ",
                     "observed block%s the biclique search extracted, and ",
                     "each was cut to exactly those %d periods. Kept %d of ",
                     "%d observations (%.1f%%). The discarded pairs and ",
                     "periods are what buys exact validity under ",
                     "missingness."),
              L0, lev_txt, L0, n_t, pd$n_mask, pd$ncell, n_cells_used,
              length(blocks), if (length(blocks) == 1L) "" else "s", L0,
              Nk, N, 100 * Nk / N)
    },
    sprintf("Block sizes (rows x cols): %s.",
            paste(sprintf("%dx%d", sizes[1, ], sizes[2, ]), collapse = ", ")),
    if (K_was_null && is.numeric(alpha) && length(alpha) == 1L &&
        is.finite(alpha) && alpha > 0 && alpha < 1 && 1 / (K + 1) > alpha)
      sprintf(paste0("Resolution here is set by the smallest selected ",
                     "block: its permuted side is %d, so K = %d. Raise ",
                     "`min_block` so that small blocks cannot set K, or ",
                     "%s so that more pairs clear the mask."),
              pd$min_side, K,
              if (is.null(L0)) "drop the sparsest periods"
              else "lower `L0`")
    else character(0))

  ## --- (iv) Procedure 2 with the period held fixed --------------------------
  perm_builder <- function(rep_seed) pd$perm_builder(rep_seed, K)

  res <- .ipt_engine(yk, Dk, Xk, perm_builder, K = K, n_reps = n_reps,
                     seed = seed, alpha = alpha, conf_int = conf_int,
                     beta_null = beta_null, grid = grid,
                     type = "panel (bicliques)", d_names = d_names,
                     n_clusters = c(row = pd$n_row, col = pd$n_col,
                                    time = length(S)),
                     call = cl, n_cores = n_cores, ci_agg = aggregate)
  res$note <- c(note, res$note)
  res$n_blocks <- length(blocks)
  res$cells_used <- n_cells_used
  res$cells_total <- pd$ncell
  res$periods_used <- pd$S_labels           # NULL when L0 is NULL
  res
}

#' The common period set under `L0`, chosen from the observation pattern.
#'
#' Given the cell x period incidence `inc` (TRUE where the pair observes the
#' period), returns L0 period indices jointly observed by as many pairs as the
#' greedy walk can find: start from the period observed by the most pairs,
#' then repeatedly add the period that keeps the most pairs observing the
#' whole set so far. Ties go to the lower period index (`which.max()` returns
#' the first maximum), so the choice is deterministic and, with periods
#' ordered by `factor()`, reproducible across calls. Only pairs with at least
#' L0 distinct periods are counted: no other pair can clear the mask, and
#' restricting to them guarantees the walk never dead-ends, because every
#' pair still alive after s steps has at least L0 - s periods outside the set.
#'
#' Maximising the joint count exactly over all L0-subsets is a maximum-
#' biclique problem in the pair x period graph, so a heuristic is used, as
#' `find_bicliques(method = "greedy")` does for the blocks. The choice only
#' affects power: it reads the observation mask and nothing else -- never
#' `y` -- so any period set it returns is a function of the mask alone, which
#' is what Assumption 4 (mask independent of the errors given the covariates)
#' needs for the cut to cost no validity.
#'
#' When the periods are within-cell ranks, period v is observed by every pair
#' with at least v observations, the counts are non-increasing in v, the walk
#' picks 1, 2, ..., L0 in turn, and the mask "observes every period in S" is
#' exactly the count mask 1{ell_ij >= L0}.
#'
#' (Moved unchanged from `R/irregular.R` in 0.4.2, where it served the retired
#' `trim = "levels"` option.)
#'
#' @param inc logical matrix, one row per cell, one column per period.
#' @param L0 the number of periods to retain.
#' @return sorted integer vector of L0 column indices of `inc`, or
#'   `integer(0)` when no cell observes L0 distinct periods.
#' @keywords internal
#' @noRd
.choose_common_levels <- function(inc, L0) {
  deep <- rowSums(inc) >= L0
  if (!any(deep) || ncol(inc) < L0) return(integer(0))
  alive <- deep
  avail <- rep(TRUE, ncol(inc))
  S <- integer(L0)
  for (s in seq_len(L0)) {
    cnt <- colSums(inc[alive, , drop = FALSE])
    cnt[!avail] <- -1
    v <- which.max(cnt)                # first maximum: the lower level on ties
    S[s] <- v
    avail[v] <- FALSE
    alive <- alive & inc[, v]
  }
  sort(S)
}

#' Steps (i)-(iii) of the incomplete-panel design: the mask, the blocks, the
#' retained observations, and the Procedure 2 builder over them.
#'
#' Everything `mwperm_panel_missing()` does to the observation pattern before
#' the engine sees a number, factored out so that the structural claims --
#' every retained pair observes every retained period, and every gather
#' vector preserves the period -- can be asserted on the front end's OWN
#' retained set and permutation builder rather than on a hand-built imitation
#' of them (see `tests/test-panel-missing.R`).
#'
#' With `L0 = NULL` the mask is "pair observed in EVERY period" and the
#' retained period set S is all of them: the 0.4.1 design, bit for bit. With
#' an integer `L0` the set S is the L0 periods jointly observed by the most
#' pairs (`.choose_common_levels()`, from the observation pattern alone), the
#' mask is "pair observes every period in S", and observations outside S are
#' dropped; `L0` equal to the number of periods reproduces `L0 = NULL`. Either
#' way the slot handed to the block builder is the period's rank within S, so
#' cell (i, j) in period s maps to cell (pi(i), sigma(j)) in the SAME period s
#' -- condition InvB, blockwise.
#'
#' @param row,col,time identifiers, one per observation.
#' @param L0 `NULL`, or the number of periods to retain (validated by the
#'   caller as a single integer >= 2; refused here if above the period count).
#' @param min_block,block_method passed to `find_bicliques()`.
#' @return a list: `idx` (sorted indices of the retained observations), `S`
#'   (dense ids of the retained periods, `seq_len(n_t)` under `L0 = NULL`),
#'   `S_labels` (their labels; `NULL` under `L0 = NULL`), `ti_k` (dense period
#'   id of each retained observation), `ncell`, `n_mask`, `blocks`,
#'   `min_side`, `n_cells_used`, `n_row`, `n_col`, `n_t`, and
#'   `perm_builder(rep_seed, K)`.
#' @keywords internal
#' @noRd
.panel_missing_design <- function(row, col, time, L0 = NULL, min_block = 3L,
                                  block_method = "greedy") {
  ri <- .dense_id(row, "row")
  ci <- .dense_id(col, "col")
  tf <- factor(time)
  ti <- .dense_id(tf, "time")
  n_row <- max(ri)
  n_col <- max(ci)
  n_t <- max(ti)

  ## One observation per (row, col, period). Without this the "same (pi, sigma)
  ## in every period" map is not well defined: two observations in one cell and
  ## period would both claim the same slot, and the gather vector would resolve
  ## the collision silently to whichever came first.
  if (anyDuplicated(cbind(ri, ci, ti)))
    stop(paste0("Each (row, col, time) cell must appear at most once. For ",
                "several observations per cell use mwperm_layout() or, if ",
                "`d` is constant within a cell, mwperm_irregular(); both ",
                "need the repeats to be exchangeable replicates."),
         call. = FALSE)

  cell <- .dense_id(interaction(ri, ci, drop = TRUE))
  ncell <- max(cell)
  first <- match(seq_len(ncell), cell)

  ## --- (i) the mask ----------------------------------------------------------
  if (is.null(L0)) {
    ## Pairs observed in EVERY period. Given uniqueness above, a cell has all
    ## n_t periods exactly when it holds n_t observations, so the count is the
    ## set condition.
    ell <- tabulate(cell, nbins = ncell)
    mask <- ell == n_t
    if (!any(mask))
      stop(sprintf(paste0("No (row, col) pair is observed in all %d periods, ",
                          "so the mask is empty (the fullest pair has %d of ",
                          "%d). Drop the sparsest periods and call again, ",
                          "use mwperm_missing() on a single period, or pass ",
                          "`L0 =` to keep the best-covered L0 periods."),
                   n_t, max(ell), n_t), call. = FALSE)
    S <- seq_len(n_t)
    S_labels <- NULL
  } else {
    if (L0 > n_t)
      stop(sprintf(paste0("`L0` = %d exceeds the number of periods (%d): ",
                          "at most %d periods can be retained. Use L0 <= ",
                          "%d, or leave L0 = NULL to keep every period."),
                   L0, n_t, n_t, n_t), call. = FALSE)
    ## Pairs observing every period of the common set S. The incidence is
    ## the observation pattern and nothing else (Assumption 4).
    inc <- matrix(FALSE, ncell, n_t)
    inc[cbind(cell, ti)] <- TRUE
    ell <- rowSums(inc)
    if (max(ell) < L0)
      stop(sprintf(paste0("No (row, col) pair is observed in L0 = %d ",
                          "periods (the fullest pair has %d of %d), so no ",
                          "common set of L0 periods exists and the mask is ",
                          "empty. Lower `L0`."),
                   L0, max(ell), n_t), call. = FALSE)
    S <- .choose_common_levels(inc, L0)
    mask <- rowSums(inc[, S, drop = FALSE]) == L0
    S_labels <- levels(tf)[S]
  }

  ## --- (ii) biclique search on the mask (Algorithm 2) -----------------------
  blocks <- find_bicliques(ri[first[mask]], ci[first[mask]],
                           min_block = min_block, method = block_method)
  if (length(blocks) == 0L)
    stop(sprintf(paste0("No fully observed block with both sides >= ",
                        "min_block = %s was found among the %d pairs ",
                        "observed in %s. Lower `min_block`, or drop the ",
                        "sparsest periods%s so that more pairs clear the ",
                        "mask."),
                 paste(min_block, collapse = "x"), sum(mask),
                 if (is.null(L0)) sprintf("all %d periods", n_t)
                 else sprintf("every retained period {%s}",
                              paste(S_labels, collapse = ", ")),
                 if (is.null(L0)) " (or pass `L0 =`)" else " (or lower `L0`)"),
         call. = FALSE)

  ## --- (iii) restrict to the cells inside the selected blocks, and to S -----
  in_block <- integer(ncell)              # block index q per cell; 0 = none
  for (q in seq_along(blocks)) {
    b <- blocks[[q]]
    hit <- mask & (ri[first] %in% b$rows) & (ci[first] %in% b$cols)
    in_block[hit] <- q
  }
  idx <- which(in_block[cell] > 0L & ti %in% S)   # retained observations
  if (!length(idx))
    stop("Internal error: the selected blocks contain no observations.",
         call. = FALSE)   # nocov
  ri_k <- ri[idx]
  ci_k <- ci[idx]
  ti_k <- ti[idx]
  Nk <- length(idx)

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
  min_side <- min(vapply(blocks,
                         function(b) min(length(b$rows), length(b$cols)),
                         integer(1)))

  ## --- (iv) Procedure 2 with the period held fixed --------------------------
  ## The slot is the period's rank within S (1..|S|): the period itself under
  ## L0 = NULL, where S is every period and match() is the identity on the
  ## dense ids.
  slot <- match(ti_k, S)
  perm_builder <- function(rep_seed, K)
    .build_obs_perms_blocks(rep_seed, K, blocks,
                            ri = ri_k, ci = ci_k, blk = blk_k,
                            lrow = lrow_k, lcol = lcol_k, permute = "both",
                            slot = slot)

  list(idx = idx, S = S, S_labels = S_labels, ti_k = ti_k,
       ncell = ncell, n_mask = sum(mask), blocks = blocks,
       min_side = min_side, n_cells_used = sum(in_block > 0L),
       n_row = n_row, n_col = n_col, n_t = n_t,
       perm_builder = perm_builder)
}
