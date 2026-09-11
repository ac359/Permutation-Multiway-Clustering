#' Invariant permutation test for irregular two-way layouts (Section 6.4)
#'
#' Finite-sample valid test of H0: beta = b for a two-way layout
#' \preformatted{  y_ijl = x_ijl' gamma + d_ijl' beta + eps_ijl}
#'
#' where cell (i, j) holds ell_ij observations indexed by l, the cell sizes
#' are unequal, and permuting *within* a cell is either invalid or powerless.
#' This is the procedure of Guo, Toulis and Wang (2026), Section 6.4. Two
#' cases motivate it, and [mwperm_layout()] handles neither:
#' - the replication index l is really *time*, so the errors are not
#'   exchangeable across l (invariance InvB, not InvA) and within-cell
#'   permutation is **invalid**;
#' - d is constant within each (i, j) cell (a dyad-level covariate), so
#'   within-cell permutation leaves the residual statistic unchanged and the
#'   test has **no power**.
#'
#' In both cases the fix is to permute the *cells*, across i and j, exactly as
#' [mwperm_dyadic()] does -- which needs equal cell sizes and a complete
#' array, neither of which an irregular layout has. Section 6.4 obtains both
#' by combining the missing-data machinery of Procedure 2 with the panel
#' construction of Case (B):
#' 1. record, for every cell, the set of `rep` levels it observes (its order
#'   of appearance when `rep` is `NULL`), and choose from that pattern alone a
#'   common set S of L0 levels -- the L0 levels jointly observed by the most
#'   cells, built up greedily with ties going to the lower level. The mask is
#'   M_ij = 1 if cell (i, j) observes every level in S, and 0 otherwise;
#' 2. run the biclique search ([find_bicliques()], Algorithm 2) on that mask
#'   to obtain disjoint fully observed blocks I_q by J_q;
#' 3. inside each selected block, keep exactly the observations at the levels
#'   in S, so every retained cell holds the same L0 levels and nothing else;
#' 4. apply Procedure 2 to what remains: a row group on I_q and a column group
#'   on J_q per block, of common order K + 1, applied identically at every
#'   level, so cell (i, j) level s maps to cell (pi(i), sigma(j)) level s.
#'
#' Holding the level fixed is the same device [mwperm_panel()] uses for time,
#' and it is what makes the test valid when the within-cell index is a period.
#'
#' @section Departure from the printed step (i): Section 6.4 as printed masks
#'   on the cell *count*, M_ij = 1 if ell_ij >= L0, and then drops ell_ij - L0
#'   observations from each retained cell at random. That leaves every cell
#'   with L0 observations but says nothing about *which* L0: after independent
#'   per-cell draws, the l-th survivor is one period in one cell and a
#'   different period in the next, and the index the permutation then holds
#'   fixed is no longer the t of condition InvB. With a covariate that varies
#'   within cells and any common period effect in the errors, that test is not
#'   valid (empirical size near 1 at nominal 0.05 in this package's own
#'   simulations). Section 6.4 is defined as Section 5 combined with Case (B),
#'   whose invariance carries the *same* t on both sides, so this function
#'   retains a common level set instead. When every cell observes the levels
#'   1..ell_ij in order -- `rep = NULL` -- the common set is `{1, ..., L0}`
#'   and the mask reduces exactly to the printed 1{ell_ij >= L0}, with the
#'   first L0 observations of each cell kept.
#'
#' @section Assumptions: This test does **not** assume within-cell
#'   exchangeability. What it needs is exchangeability across the cell indices
#'   (i, j) *within each slot*: for every slot l, the errors in a block
#'   satisfy
#'
#' (eps_ijl for i in I_q, j in J_q) has the same distribution, given X and D,
#' as (eps_[pi(i)][sigma(j)]l for i in I_q, j in J_q)
#'
#' with the *same* (pi, sigma) used in every slot -- invariance (InvB) of
#' Section 6.2, applied blockwise. That is required together with Assumption 4
#' on the mask: the missingness pattern (here, which cells observe the common
#' level set S) is independent of the errors given the covariates. Under
#' those two conditions the p-value is exact in finite samples by the
#' argument of Theorem 4. Arbitrary slot effects shared across cells -- an
#' arbitrary common time trend -- are permitted, which is the whole point:
#' every retained observation keeps its `rep` level under every permutation,
#' so a level effect zeta_l is carried along unchanged, exactly as the period
#' effect is in [mwperm_panel()].
#'
#' @section Choosing L0: L0 trades cells against within-cell depth: a small L0
#'   keeps more cells in the mask (so larger blocks, larger K, finer p-value
#'   resolution) but throws away more observations per cell; a large L0 keeps
#'   deeper cells but fewer of them. The paper recommends tuning it by grid
#'   search to minimise the loss of observations. There is no free lunch in
#'   choosing it from the data on the *outcome*; select it from the
#'   observation pattern alone, which is ancillary under Assumption 4. The
#'   common level set S is chosen the same way, from the pattern of which cell
#'   observes which level and nothing else; which levels it picks affects
#'   power, never validity. To dictate the levels yourself, subset the data to
#'   the levels you want and set `L0` to their number: S is then exactly that
#'   set. The `note` field and the `rep_levels` field of the fitted object
#'   report which levels, cells and observations survived.
#'
#' @section Reproducibility: The retained observations are a deterministic
#'   function of the observation pattern and `L0`; nothing is subsampled at
#'   random. The only randomness is the permutation group, so `seed` governs
#'   the fit exactly as it does in [mwperm_missing()].
#'
#' @inheritParams mwperm_dyadic
#' @param d Numeric vector or matrix of the covariate(s) of interest. Unlike
#'   [mwperm_layout()], `d` *may* be constant within cells; that is one of the
#'   cases this design exists for, and no warning is issued.
#' @param x Optional numeric matrix or data frame of nuisance covariates; an
#'   intercept is always added internally. May be `NULL`.
#' @param row,col Cell identifiers along the two layout dimensions.
#' @param rep Optional within-cell index (the period or wave identifier). It
#'   is the index the permutation holds fixed: cell (i, j) level s maps to
#'   cell (pi(i), sigma(j)) level s, and a cell is retained only if it
#'   observes every level in the common set S, which is what makes level s
#'   the same period in every cell. It must be unique within a cell. When
#'   `NULL`, order of appearance within the cell is used, and the mask reduces
#'   to the printed 1{ell_ij >= L0}.
#' @param L0 Integer, at least 2: the number of `rep` levels retained per
#'   cell. A common set S of L0 levels is chosen from the observation pattern,
#'   cells that do not observe every level in S are masked out, and every
#'   retained cell is cut to exactly those L0 observations. Required.
#' @param K Number of non-identity permutations; defaults to the smallest
#'   block side over the selected blocks -- that is, the smallest of |I_q| and
#'   |J_q| over all q -- minus one, capped at 199. Must satisfy `K + 1 <=`
#'   that smallest block side.
#' @param min_block Minimum block side(s) for the biclique search; see
#'   [find_bicliques()].
#' @param block_method `"greedy"` (default) or `"exact"`; see
#'   [find_bicliques()].
#'
#' @return An object of class `"mwperm"`, with the extra fields `n_blocks`,
#'   `cells_used`, `cells_total`, `L0` and `rep_levels` (the labels of the L0
#'   retained levels). `estimate`/`se_naive` are the OLS
#'   estimate and naive SE on the retained data, `conf_int`/`conf_set` (or
#'   `conf_region`/`conf_box` for several coefficients) the IPT inverted-test
#'   confidence set, and `pvalue` the IPT permutation p-value; see
#'   [mwperm_dyadic()] for the field provenance in full.
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data, Section 6.4 and Procedure 2.
#'   arXiv:2601.08610.
#' @seealso [mwperm_layout()] (Section 6.3, within-cell permutation),
#'   [mwperm_missing()] (Procedure 2 with one observation per cell),
#'   [mwperm_panel()] (Section 6.2, complete balanced panels).
#' @examples
#' ## 8 x 8 cells, unequal cell sizes, treatment CONSTANT within each cell --
#' ## the case where the within-cell test of mwperm_layout() has no power.
#' set.seed(11)
#' dat <- do.call(rbind, lapply(1:8, function(i)
#'   do.call(rbind, lapply(1:8, function(j) {
#'     L <- sample(c(0L, 2L, 4L, 6L, 9L), 1L)
#'     if (L == 0L) return(NULL)
#'     data.frame(i = i, j = j, l = seq_len(L), d = rnorm(1), eta = rnorm(1))
#'   }))))
#' dat$y <- 0.5 * dat$d + dat$eta + rnorm(nrow(dat))
#' with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
#'                            L0 = 4L, min_block = 2L, conf_int = FALSE,
#'                            seed = 1))
#' @export
mwperm_irregular <- function(y, d, x = NULL, row, col, rep = NULL, L0,
                             K = NULL, min_block = 3L,
                             block_method = c("greedy", "exact"),
                             alpha = 0.05, beta_null = 0, conf_int = TRUE,
                             n_reps = 10L, seed = NULL, grid = NULL,
                             aggregate = c("median", "median2"),
                             n_cores = 1L) {
  cl <- match.call()
  block_method <- match.arg(block_method)
  aggregate <- match.arg(aggregate)
  y <- .check_y(y)
  N <- length(y)
  D <- as.matrix(d)
  d_names <- .coef_names(D, deparse(substitute(d)))
  X <- .make_X(x, N)
  .check_lengths(N, list(row = row, col = col,
                         rep = if (is.null(rep)) seq_len(N) else rep))
  if (!is.null(rep) && anyNA(rep))
    stop(paste0("`rep` contains missing values (NA); the within-cell ",
                "identifier must be complete."), call. = FALSE)
  ## Validate the FULL supplied data before the mask and the biclique step
  ## discard anything: the observed cells are the user's data contract.
  .check_finite(list(y = y, d = D, x = X))

  if (missing(L0))
    stop(paste0("`L0` is required: it is the number of `rep` levels retained ",
                "per cell, which sets the common level set S and the mask ",
                "M_ij = 1{cell (i, j) observes every level in S} (Section ",
                "6.4). Pick it from the observation pattern -- see the ",
                "`Choosing L0` section of ?mwperm_irregular."), call. = FALSE)
  L0 <- suppressWarnings(as.integer(L0))
  if (length(L0) != 1L || is.na(L0) || L0 < 2L)
    stop("`L0` must be a single integer >= 2.", call. = FALSE)

  ## --- (i)-(iii) the aligned mask, the blocks and the deterministic cut ------
  pd <- .irregular_design(row, col, rep, L0, min_block = min_block,
                          block_method = block_method)
  idx <- pd$idx                              # retained observations, sorted
  yk <- y[idx]
  Dk <- D[idx, , drop = FALSE]
  Xk <- X[idx, , drop = FALSE]
  Nk <- length(idx)
  blocks <- pd$blocks
  n_cells_used <- pd$n_cells_used

  ## --- group order: the smallest permuted block side over the blocks ---------
  K_was_null <- is.null(K)
  K <- .default_K(K, pd$min_side)

  sizes <- vapply(blocks, function(b) c(length(b$rows), length(b$cols)),
                  integer(2))
  lev_txt <- if (L0 <= 10L) paste(pd$S_labels, collapse = ", ")
             else sprintf("%s, ..., %s",
                          paste(pd$S_labels[1:3], collapse = ", "),
                          pd$S_labels[L0])
  note <- c(
    sprintf(paste0("Section 6.4 (irregular design), L0 = %d: the retained ",
                   "levels are {%s}%s, the %d levels jointly observed by the ",
                   "most cells; %d of %d cells observe all of them, %d of ",
                   "those lie in the %d fully observed block%s the biclique ",
                   "search extracted, and each was cut to exactly those %d ",
                   "observations. Kept %d of %d observations (%.1f%%). The ",
                   "discarded observations are what buys exact validity ",
                   "under an unequal design."),
            L0, lev_txt,
            if (is.null(rep)) " (order of appearance within the cell)" else "",
            L0, pd$n_mask, pd$ncell, n_cells_used, length(blocks),
            if (length(blocks) == 1L) "" else "s", L0, Nk, N, 100 * Nk / N),
    sprintf("Block sizes (rows x cols): %s.",
            paste(sprintf("%dx%d", sizes[1, ], sizes[2, ]), collapse = ", ")),
    if (K_was_null && is.numeric(alpha) && length(alpha) == 1L &&
        is.finite(alpha) && alpha > 0 && alpha < 1 && 1 / (K + 1) > alpha)
      sprintf(paste0(
        "Resolution here is set by the smallest selected block: its permuted ",
        "side is %d, so K = %d. Raise `min_block` so that small blocks ",
        "cannot set K, or retune L0 -- a larger L0 gives deeper cells but ",
        "fewer of them, a smaller one admits more cells and can support ",
        "larger blocks."),
        pd$min_side, K)
    else character(0))

  ## --- (iv) Procedure 2 on the stacked retained data, level held fixed ------
  perm_builder <- function(rep_seed) pd$perm_builder(rep_seed, K)

  res <- .ipt_engine(yk, Dk, Xk, perm_builder, K = K, n_reps = n_reps,
                     seed = seed, alpha = alpha, conf_int = conf_int,
                     beta_null = beta_null, grid = grid,
                     type = "irregular (Section 6.4)", d_names = d_names,
                     n_clusters = c(row = pd$n_row, col = pd$n_col),
                     call = cl, n_cores = n_cores, ci_agg = aggregate)
  res$note <- c(note, res$note)
  res$n_blocks <- length(blocks)
  res$cells_used <- n_cells_used
  res$cells_total <- pd$ncell
  res$L0 <- L0
  res$rep_levels <- pd$S_labels
  res
}

#' The common level set of Section 6.4, chosen from the observation pattern.
#'
#' Given the cell x level incidence `inc` (TRUE where the cell observes the
#' level), returns L0 level indices jointly observed by as many cells as the
#' greedy walk can find: start from the level observed by the most cells,
#' then repeatedly add the level that keeps the most cells observing the whole
#' set so far. Ties go to the lower level index (`which.max()` returns the
#' first maximum), so the choice is deterministic and, with levels ordered by
#' `factor()`, reproducible across calls. Only cells with at least L0 distinct
#' levels are counted: no other cell can clear the mask, and restricting to
#' them guarantees the walk never dead-ends, because every cell still alive
#' after s steps has at least L0 - s levels outside the set.
#'
#' Maximising the joint count exactly over all L0-subsets is a maximum-
#' biclique problem in the cell x level graph, so a heuristic is used, as
#' `find_bicliques(method = "greedy")` does for the blocks. The choice only
#' affects power: any level set chosen from the observation pattern alone is
#' ancillary under Assumption 4.
#'
#' When the levels are within-cell ranks (`rep = NULL`), level v is observed
#' by every cell with ell_ij >= v, so the counts are non-increasing in v, the
#' walk picks 1, 2, ..., L0 in turn, and the mask "observes every level in S"
#' is exactly GTW's printed 1{ell_ij >= L0}. `tests/test-irregular.R` pins
#' that degeneracy.
#'
#' @param inc logical matrix, one row per cell, one column per level.
#' @param L0 the number of levels to retain.
#' @return sorted integer vector of L0 column indices of `inc`, or
#'   `integer(0)` when no cell observes L0 distinct levels.
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

#' Steps (i)-(iii) of Section 6.4: the aligned mask, the blocks, the cut.
#'
#' Everything `mwperm_irregular()` does to the observation pattern before the
#' engine sees a number, factored out so that the structural claim -- every
#' gather vector preserves the `rep` level -- can be asserted on the front
#' end's OWN retained set and permutation builder rather than on a hand-built
#' imitation of them (see `tests/test-irregular.R`).
#'
#' @param row,col cell identifiers, one per observation.
#' @param rep within-cell level identifier, or `NULL` for order of appearance.
#' @param L0 number of levels retained per cell.
#' @param min_block,block_method passed to `find_bicliques()`.
#' @return a list: `idx` (sorted indices of the retained observations),
#'   `slot` (their level, re-indexed 1..L0 in level order, the coordinate the
#'   permutation holds fixed), `S` and `S_labels` (the retained level indices
#'   and labels), `cell`, `ncell`, `n_mask`, `in_block` (block of each cell,
#'   0 = discarded), `blocks` (those that retained cells), `min_side`,
#'   `n_cells_used`, `n_row`, `n_col`, and `perm_builder(rep_seed, K)`, the
#'   Procedure 2 builder over the retained observations.
#' @keywords internal
#' @noRd
.irregular_design <- function(row, col, rep = NULL, L0, min_block = 3L,
                              block_method = "greedy") {
  ri <- .dense_id(row, "row")
  ci <- .dense_id(col, "col")
  cell <- .dense_id(interaction(ri, ci, drop = TRUE))   # dense occupied-cell id
  ncell <- max(cell)

  ## --- (i) the level of every observation, and the common set S -------------
  ## The level is the coordinate the permutation will hold fixed, so it has to
  ## mean the same thing in every cell: the user's `rep` where given, the order
  ## of appearance otherwise. Two observations of one cell at one level would
  ## both claim the same (row, col, level) key and the gather vector would
  ## resolve the collision silently to the first, so that is refused up front,
  ## exactly as mwperm_panel_missing() refuses a repeated (row, col, time).
  if (is.null(rep)) {
    lev <- .within_cell_slot(cell, NULL, ncell)
    S_labels_all <- as.character(seq_len(max(lev)))
  } else {
    f <- factor(rep)
    lev <- .dense_id(f, "rep")
    S_labels_all <- levels(f)
    dup <- anyDuplicated(cbind(cell, lev))
    if (dup > 0L)
      stop(sprintf(paste0("Cell (row = %s, col = %s) holds `rep` level %s ",
                          "more than once (observation %d). `rep` is the ",
                          "index the permutation holds fixed, so it must be ",
                          "unique within a cell: if the repeats are ",
                          "exchangeable replicates rather than periods, use ",
                          "mwperm_layout(), or drop `rep` to index them by ",
                          "order of appearance."),
                   row[dup], col[dup], S_labels_all[lev[dup]], dup),
           call. = FALSE)
  }
  n_lev <- max(lev)
  inc <- matrix(FALSE, ncell, n_lev)
  inc[cbind(cell, lev)] <- TRUE                         # cell observes level
  ell <- rowSums(inc)                                   # distinct levels/cell
  if (max(ell) < L0)
    stop(sprintf(paste0("No cell has at least L0 = %d distinct `rep` levels ",
                        "(the deepest has %d), so no common set of L0 levels ",
                        "exists and the mask is empty. Lower `L0`."),
                 L0, max(ell)), call. = FALSE)
  S <- .choose_common_levels(inc, L0)
  mask <- rowSums(inc[, S, drop = FALSE]) == L0         # observes all of S
  first <- match(seq_len(ncell), cell)                  # a row of each cell

  ## --- (ii) biclique search on the mask (Algorithm 2) ------------------------
  blocks <- find_bicliques(ri[first[mask]], ci[first[mask]],
                           min_block = min_block, method = block_method)
  if (length(blocks) == 0L)
    stop(sprintf(paste0("No fully observed block with both sides >= ",
                        "min_block = %s was found among the %d cells that ",
                        "observe every retained level {%s}. Lower ",
                        "`min_block`, lower `L0` (a smaller level set is ",
                        "observed by more cells), or use mwperm_layout() if ",
                        "within-cell permutation is valid for these data."),
                 paste(min_block, collapse = "x"), sum(mask),
                 paste(S_labels_all[S], collapse = ", ")),
         call. = FALSE)

  ## --- restrict to cells inside the selected blocks --------------------------
  ## Cells outside every block are dropped whole; that discard is the price of
  ## validity under Procedure 2, exactly as in mwperm_missing().
  in_block <- integer(ncell)          # block index q of each cell; 0 = none
  for (q in seq_along(blocks)) {
    b <- blocks[[q]]
    hit <- mask & (ri[first] %in% b$rows) & (ci[first] %in% b$cols)
    in_block[hit] <- q
  }

  ## --- (iii) keep exactly the observations at the levels in S ----------------
  ## Deterministic: a retained cell observes every level in S (the mask) and
  ## each at most once (the check above), so it keeps exactly L0 observations,
  ## and level s is the same level in every cell.
  idx <- which(in_block[cell] > 0L & lev %in% S)
  if (!length(idx))
    stop("Internal error: the selected blocks contain no observations.",
         call. = FALSE)   # nocov
  slot <- match(lev[idx], S)                 # 1..L0, in level order
  ri_k <- ri[idx]
  ci_k <- ci[idx]
  cell_k <- cell[idx]
  n_cells_used <- length(unique(cell_k))

  ## Block index and block-local row/column position of every retained
  ## observation, for the Procedure 2 permutation builder.
  blk_k <- in_block[cell_k]
  lrow_k <- integer(length(idx))
  lcol_k <- integer(length(idx))
  for (q in seq_along(blocks)) {
    s <- which(blk_k == q)
    if (!length(s)) next
    lrow_k[s] <- match(ri_k[s], blocks[[q]]$rows)
    lcol_k[s] <- match(ci_k[s], blocks[[q]]$cols)
  }
  ## Blocks that actually kept cells (defensive: under the mask every block
  ## cell clears S, so this is the identity, but the group order below must
  ## never be set by a block with nothing in it).
  used_q <- sort(unique(blk_k))
  blocks <- blocks[used_q]
  blk_k <- match(blk_k, used_q)
  min_side <- min(vapply(blocks,
                         function(b) min(length(b$rows), length(b$cols)),
                         integer(1)))

  perm_builder <- function(rep_seed, K)
    .build_obs_perms_blocks(rep_seed, K, blocks,
                            ri = ri_k, ci = ci_k, blk = blk_k,
                            lrow = lrow_k, lcol = lcol_k, permute = "both",
                            slot = slot)

  list(idx = idx, slot = slot, S = S, S_labels = S_labels_all[S],
       cell = cell, ncell = ncell, ell = ell, n_mask = sum(mask),
       in_block = in_block, blocks = blocks, blk = blk_k,
       min_side = min_side, n_cells_used = n_cells_used,
       n_row = max(ri), n_col = max(ci), perm_builder = perm_builder)
}
