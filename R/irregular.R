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
#' construction of Case (B). As printed, and as the paper applies it in its
#' Appendix B (`trim = "random"`, the default):
#' 1. form the cell sizes ell_ij and, for a threshold L0, the mask M_ij = 1 if
#'   ell_ij >= L0, and 0 otherwise;
#' 2. run the biclique search ([find_bicliques()], Algorithm 2) on that mask
#'   to obtain disjoint fully observed blocks I_q by J_q;
#' 3. inside each selected block, drop ell_ij - L0 observations from each
#'   cell **uniformly at random**, leaving exactly L0 everywhere;
#' 4. apply Procedure 2 to what remains: a row group on I_q and a column group
#'   on J_q per block, of common order K + 1, applied identically across the
#'   L0 within-cell positions, so the l-th remaining observation of cell
#'   (i, j) maps to the l-th remaining observation of cell (pi(i), sigma(j)).
#'
#' Step 3 is random, so -- as the paper notes at the end of Section 6.4 and
#' does in Appendix B -- the whole test is repeated and the median p-value
#' reported: every one of the `n_reps` repetitions here draws its **own**
#' subsample from its own seed, together with its own permutation group, and
#' the reported p-value and confidence set aggregate over both (see
#' *Aggregation over repetitions* in [mwperm_dyadic()]).
#'
#' @section When the repeats are periods (`trim = "levels"`): The random trim
#'   leaves each cell with L0 observations but says nothing about *which*
#'   L0. When the within-cell index is a period with a common effect zeta_l
#'   in the errors, condition InvB carries the *same* l on both sides of the
#'   invariance, so the retained observations must be the same L0 periods in
#'   every cell; after independent per-cell draws the l-th survivor is one
#'   period in one cell and another period in the next, and the position the
#'   permutation holds fixed is no longer the l of InvB (with a covariate
#'   that varies within cells and a strong common period effect, that test
#'   rejected a true null almost always in this package's simulations). For
#'   that case set `trim = "levels"`: a common set S of L0 `rep` levels is
#'   chosen from the observation pattern alone (the L0 levels jointly
#'   observed by the most cells, greedily, ties to the lower level), the mask
#'   becomes M_ij = 1 if cell (i, j) observes every level in S, and exactly
#'   the observations at those levels are kept -- a deterministic cut, with
#'   level s the same period in every cell. When every cell observes the
#'   levels 1..ell_ij in order (`rep = NULL`) the common set is
#'   `{1, ..., L0}` and the mask is again the printed 1{ell_ij >= L0}, with
#'   the first L0 observations of each cell kept.
#'
#' @section Assumptions: This test does **not** assume within-cell
#'   exchangeability. What it needs is exchangeability of the retained error
#'   array across the cell indices (i, j) within each block, position by
#'   position: with the *same* (pi, sigma) at every within-cell position l,
#'
#' (eps_ijl for i in I_q, j in J_q) has the same distribution, given X and D,
#' as (eps_[pi(i)][sigma(j)]l for i in I_q, j in J_q),
#'
#' together with Assumption 4 on the mask: which cells clear L0 (and, under
#' `trim = "levels"`, which observe the common level set) is independent of
#' the errors given the covariates. Under those two conditions the p-value is
#' exact in finite samples by the argument of Theorem 4. Under
#' `trim = "random"` the position l of a retained observation is its rank
#' among the survivors, so the condition holds when the repeats inside a cell
#' are exchangeable replicates -- individuals sampled within a cell, as in
#' the paper's Appendix B -- with any cell-level dependence of the
#' random-effects kind (eps_ijl = eta_i + xi_j + u_ijl). It does **not** hold
#' under a common effect indexed by l, an arbitrary period effect zeta_l,
#' unless the retained l are aligned across cells; that is what
#' `trim = "levels"` guarantees, and under it zeta_l is carried along
#' unchanged exactly as the period effect is in [mwperm_panel()].
#'
#' @section Choosing L0: L0 trades cells against within-cell depth: a small L0
#'   keeps more cells in the mask (so larger blocks, larger K, finer p-value
#'   resolution) but throws away more observations per cell; a large L0 keeps
#'   deeper cells but fewer of them. The paper recommends tuning it by grid
#'   search to minimise the loss of observations. There is no free lunch in
#'   choosing it from the data on the *outcome*; select it from the cell
#'   sizes alone, which are ancillary under Assumption 4. The `note` field of
#'   the fitted object reports how many cells and observations survived, and
#'   under `trim = "levels"` the `rep_levels` field reports which levels.
#'
#' @section Reproducibility: Under `trim = "random"` the subsample of each
#'   repetition is drawn from that repetition's seed (rep r uses
#'   `seed + r - 1`), with the usual RNG-state hygiene: a seeded call leaves
#'   the caller's random stream exactly as it found it, and the rep-parallel
#'   path (`n_cores > 1` with a `seed`) reproduces the serial result bit for
#'   bit. Under `trim = "levels"` the retained observations are a
#'   deterministic function of the observation pattern and `L0`, and the only
#'   randomness is the permutation group.
#'
#' @inheritParams mwperm_dyadic
#' @param d Numeric vector or matrix of the covariate(s) of interest. Unlike
#'   [mwperm_layout()], `d` *may* be constant within cells; that is one of the
#'   cases this design exists for, and no warning is issued.
#' @param x Optional numeric matrix or data frame of nuisance covariates; an
#'   intercept is always added internally. May be `NULL`.
#' @param row,col Cell identifiers along the two layout dimensions.
#' @param rep Optional within-cell index (the replication or period
#'   identifier). Under `trim = "random"` it only orders the observations
#'   inside a cell (the survivors are ranked by it; order of appearance when
#'   `NULL`), so the fit does not depend on it beyond that order. Under
#'   `trim = "levels"` it is the level the permutation holds fixed -- cell
#'   (i, j) level s maps to cell (pi(i), sigma(j)) level s -- and must be
#'   unique within a cell; when `NULL`, order of appearance is used.
#' @param L0 Integer, at least 2: the number of observations retained per
#'   cell. Cells with fewer than `L0` observations (under `trim = "levels"`,
#'   cells that do not observe every level of the common set) are masked
#'   out, and every retained cell is reduced to exactly `L0`. Required.
#' @param trim How each retained cell is reduced to exactly `L0`
#'   observations. `"random"` (the default) is step (i) of Section 6.4 as
#'   printed: drop `ell_ij - L0` observations at random, independently in
#'   every cell, redrawn in every repetition. `"levels"` keeps the same L0
#'   `rep` levels in every cell, which condition InvB requires when the
#'   repeats are periods with a common effect; see *When the repeats are
#'   periods*.
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
#'   `cells_used`, `cells_total`, `L0`, `trim` and `rep_levels` (the labels
#'   of the L0 retained levels under `trim = "levels"`; `NULL` otherwise).
#'   `n_obs` is the number of observations each repetition is computed on,
#'   `cells_used * L0`. `estimate`/`se_naive` are the OLS estimate and naive
#'   SE on the retained data -- under `trim = "random"`, on *all*
#'   observations of the retained cells, before the per-repetition trim, so
#'   that the point estimate does not depend on a subsampling draw;
#'   `conf_int`/`conf_set` (or `conf_region`/`conf_box` for several
#'   coefficients) the IPT inverted-test confidence set, and `pvalue` the
#'   IPT permutation p-value; see [mwperm_dyadic()] for the field provenance
#'   in full.
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data, Section 6.4, Appendix B and
#'   Procedure 2. arXiv:2601.08610.
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
#' ## the repeats are periods with a common effect: keep the SAME periods in
#' ## every cell instead of a random subset
#' with(dat, mwperm_irregular(y = y, d = d, row = i, col = j, rep = l,
#'                            L0 = 4L, trim = "levels", min_block = 2L,
#'                            conf_int = FALSE, seed = 1))
#' @export
mwperm_irregular <- function(y, d, x = NULL, row, col, rep = NULL, L0,
                             trim = c("random", "levels"),
                             K = NULL, min_block = 3L,
                             block_method = c("greedy", "exact"),
                             alpha = 0.05, beta_null = 0, conf_int = TRUE,
                             n_reps = 10L, seed = NULL, grid = NULL,
                             aggregate = c("median", "median2"),
                             n_cores = 1L) {
  cl <- match.call()
  trim <- match.arg(trim)
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
    stop(paste0("`L0` is required: it is the cell-size threshold that ",
                "defines the mask M_ij = 1{ell_ij >= L0} and the number of ",
                "observations every retained cell is reduced to (Section ",
                "6.4). Pick it from the cell sizes -- see the `Choosing L0` ",
                "section of ?mwperm_irregular."), call. = FALSE)
  L0 <- suppressWarnings(as.integer(L0))
  if (length(L0) != 1L || is.na(L0) || L0 < 2L)
    stop("`L0` must be a single integer >= 2.", call. = FALSE)

  ## --- (i)-(iii) the mask, the blocks and the cut ----------------------------
  pd <- .irregular_design(row, col, rep, L0, trim = trim,
                          min_block = min_block, block_method = block_method)
  idx <- pd$idx                              # observations of the retained cells
  yk <- y[idx]
  Dk <- D[idx, , drop = FALSE]
  Xk <- X[idx, , drop = FALSE]
  blocks <- pd$blocks
  n_cells_used <- pd$n_cells_used
  Nk <- n_cells_used * L0                    # observations per repetition
  ## The engine checks N > 2p on the data it is handed, which under the
  ## random trim is every observation of the retained cells; the projection
  ## in each repetition sees only cells_used x L0 of them, so that is the
  ## count that must clear 2p.
  p <- ncol(Xk)
  if (Nk <= 2L * p)
    stop(sprintf(paste0("Need N > 2p for the projection to exist: each ",
                        "repetition uses N = %d retained observations ",
                        "(%d cells x L0 = %d) against p = %d nuisance ",
                        "columns (incl. intercept). Drop covariates, lower ",
                        "`min_block`, or retune L0."),
                 Nk, n_cells_used, L0, p), call. = FALSE)

  ## --- group order: the smallest permuted block side over the blocks ---------
  K_was_null <- is.null(K)
  K <- .default_K(K, pd$min_side)

  sizes <- vapply(blocks, function(b) c(length(b$rows), length(b$cols)),
                  integer(2))
  note <- c(
    if (trim == "random")
      sprintf(paste0("Section 6.4 (irregular design), L0 = %d: %d of %d ",
                     "cells have at least L0 observations, %d of those lie ",
                     "in the %d fully observed block%s the biclique search ",
                     "extracted, and in every repetition each of them is cut ",
                     "to exactly %d observations drawn at random (a fresh ",
                     "draw per repetition; the reported p-value is the ",
                     "median over %d). Each repetition uses %d of the %d ",
                     "observations (%.1f%%); the OLS estimate uses all %d ",
                     "observations of the retained cells. The discarded ",
                     "observations are what buys exact validity under an ",
                     "unequal design."),
              L0, pd$n_mask, pd$ncell, n_cells_used, length(blocks),
              if (length(blocks) == 1L) "" else "s", L0, n_reps, Nk, N,
              100 * Nk / N, length(idx))
    else {
      lev_txt <- if (L0 <= 10L) paste(pd$S_labels, collapse = ", ")
                 else sprintf("%s, ..., %s",
                              paste(pd$S_labels[1:3], collapse = ", "),
                              pd$S_labels[L0])
      sprintf(paste0("Section 6.4 (irregular design), L0 = %d, trim = ",
                     "\"levels\": the retained levels are {%s}%s, the %d ",
                     "levels jointly observed by the most cells; %d of %d ",
                     "cells observe all of them, %d of those lie in the %d ",
                     "fully observed block%s the biclique search extracted, ",
                     "and each was cut to exactly those %d observations. ",
                     "Kept %d of %d observations (%.1f%%). The discarded ",
                     "observations are what buys exact validity under an ",
                     "unequal design."),
              L0, lev_txt,
              if (is.null(rep)) " (order of appearance within the cell)"
              else "",
              L0, pd$n_mask, pd$ncell, n_cells_used, length(blocks),
              if (length(blocks) == 1L) "" else "s", L0, Nk, N, 100 * Nk / N)
    },
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

  ## --- (iv) Procedure 2 on the stacked retained data, position held fixed --
  perm_builder <- function(rep_seed) pd$perm_builder(rep_seed, K)

  res <- .ipt_engine(yk, Dk, Xk, perm_builder, K = K, n_reps = n_reps,
                     seed = seed, alpha = alpha, conf_int = conf_int,
                     beta_null = beta_null, grid = grid,
                     type = "irregular (Section 6.4)", d_names = d_names,
                     n_clusters = c(row = pd$n_row, col = pd$n_col),
                     call = cl, n_cores = n_cores, ci_agg = aggregate)
  res$note <- c(note, res$note)
  res$n_obs <- Nk                            # per repetition (= N under "levels")
  res$n_blocks <- length(blocks)
  res$cells_used <- n_cells_used
  res$cells_total <- pd$ncell
  res$L0 <- L0
  res$trim <- trim
  res$rep_levels <- pd$S_labels              # NULL under the random trim
  res
}

#' The common level set of Section 6.4 under `trim = "levels"`, chosen from
#' the observation pattern.
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

#' Steps (i)-(iii) of Section 6.4: the mask, the blocks, the cut.
#'
#' Everything `mwperm_irregular()` does to the observation pattern before the
#' engine sees a number, factored out so that the structural claims -- under
#' `trim = "random"` every repetition keeps exactly L0 observations per
#' retained cell and every gather vector preserves the within-cell position;
#' under `trim = "levels"` every gather vector preserves the `rep` level --
#' can be asserted on the front end's OWN retained set and permutation
#' builder rather than on a hand-built imitation of them (see
#' `tests/test-irregular.R`).
#'
#' Under `trim = "random"` the mask is the printed count mask
#' 1{ell_ij >= L0}, `idx` is EVERY observation of the retained cells, and the
#' cut happens inside `perm_builder(rep_seed, K)`: it draws that repetition's
#' subsample (`.downsample_to_L0()`, seeded from the rep seed at sub-seed
#' offset 1 -- the block builder uses offsets 4q - 1 and 4q, so nothing
#' collides), ranks the survivors inside each cell by `rep` (order of
#' appearance when `NULL`), builds the block-diagonal Procedure 2 group over
#' the survivors with that rank as the slot held fixed, and returns it with
#' the attribute `"rows"` naming the survivors within `idx`, which is how the
#' engine learns to compute that repetition on the subsample. Under
#' `trim = "levels"` the cut is deterministic, `idx` is already the retained
#' set, and the builder returns a plain group over it.
#'
#' @param row,col cell identifiers, one per observation.
#' @param rep within-cell level identifier, or `NULL` for order of appearance.
#' @param L0 number of observations retained per cell.
#' @param trim `"random"` or `"levels"`; see [mwperm_irregular()].
#' @param min_block,block_method passed to `find_bicliques()`.
#' @return a list: `idx` (sorted indices of the observations handed to the
#'   engine), `S` and `S_labels` (the retained level indices and labels under
#'   `"levels"`; `integer(0)` and `NULL` under `"random"`), `cell`, `ncell`,
#'   `ell` (observations per cell, or distinct levels per cell under
#'   `"levels"`), `n_mask`, `in_block` (block of each cell, 0 = discarded),
#'   `blocks` (those that retained cells), `min_side`, `n_cells_used`,
#'   `n_row`, `n_col`, `trim`, and `perm_builder(rep_seed, K)`, the
#'   Procedure 2 builder described above.
#' @keywords internal
#' @noRd
.irregular_design <- function(row, col, rep = NULL, L0,
                              trim = c("random", "levels"),
                              min_block = 3L, block_method = "greedy") {
  trim <- match.arg(trim)
  ri <- .dense_id(row, "row")
  ci <- .dense_id(col, "col")
  N <- length(ri)
  cell <- .dense_id(interaction(ri, ci, drop = TRUE))   # dense occupied-cell id
  ncell <- max(cell)

  ## --- (i) the mask ----------------------------------------------------------
  if (trim == "random") {
    ## The printed mask: cell (i, j) clears the threshold when it holds at
    ## least L0 observations, whatever their labels. `rep` only orders the
    ## observations inside a cell (the survivors are ranked by it).
    ell <- tabulate(cell, nbins = ncell)                # observations per cell
    if (max(ell) < L0)
      stop(sprintf(paste0("No cell has at least L0 = %d observations (the ",
                          "largest cell has %d), so the mask M_ij = 1{ell_ij ",
                          ">= L0} is empty. Lower `L0`."),
                   L0, max(ell)), call. = FALSE)
    mask <- ell >= L0
    ord_key <- if (is.null(rep)) seq_len(N) else as.numeric(factor(rep))
    lev <- NULL
    S <- integer(0)
    S_labels_all <- NULL
  } else {
    ## The level of every observation, and the common set S. The level is the
    ## coordinate the permutation will hold fixed, so it has to mean the same
    ## thing in every cell: the user's `rep` where given, the order of
    ## appearance otherwise. Two observations of one cell at one level would
    ## both claim the same (row, col, level) key and the gather vector would
    ## resolve the collision silently to the first, so that is refused up
    ## front, exactly as mwperm_panel_missing() refuses a repeated (row, col,
    ## time).
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
                            "more than once (observation %d). Under trim = ",
                            "\"levels\" `rep` is the index the permutation ",
                            "holds fixed, so it must be unique within a ",
                            "cell: if the repeats are exchangeable ",
                            "replicates rather than periods, use the ",
                            "default trim = \"random\" (or mwperm_layout()), ",
                            "or drop `rep` to index them by order of ",
                            "appearance."),
                     row[dup], col[dup], S_labels_all[lev[dup]], dup),
             call. = FALSE)
    }
    n_lev <- max(lev)
    inc <- matrix(FALSE, ncell, n_lev)
    inc[cbind(cell, lev)] <- TRUE                       # cell observes level
    ell <- rowSums(inc)                                 # distinct levels/cell
    if (max(ell) < L0)
      stop(sprintf(paste0("No cell has at least L0 = %d distinct `rep` ",
                          "levels (the deepest has %d), so no common set of ",
                          "L0 levels exists and the mask is empty. Lower ",
                          "`L0`."),
                   L0, max(ell)), call. = FALSE)
    S <- .choose_common_levels(inc, L0)
    mask <- rowSums(inc[, S, drop = FALSE]) == L0       # observes all of S
    ord_key <- NULL
  }
  first <- match(seq_len(ncell), cell)                  # a row of each cell

  ## --- (ii) biclique search on the mask (Algorithm 2) ------------------------
  blocks <- find_bicliques(ri[first[mask]], ci[first[mask]],
                           min_block = min_block, method = block_method)
  if (length(blocks) == 0L)
    stop(sprintf(paste0("No fully observed block with both sides >= ",
                        "min_block = %s was found among the %d cells that %s. ",
                        "Lower `min_block`, lower `L0` (%s), or use ",
                        "mwperm_layout() if within-cell permutation is valid ",
                        "for these data."),
                 paste(min_block, collapse = "x"), sum(mask),
                 if (trim == "random") sprintf("clear L0 = %d", L0)
                 else sprintf("observe every retained level {%s}",
                              paste(S_labels_all[S], collapse = ", ")),
                 if (trim == "random") "more cells clear the threshold"
                 else "a smaller level set is observed by more cells"),
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

  ## --- (iii) the cut ---------------------------------------------------------
  ## "random": every observation of a retained cell is handed to the engine
  ## and the per-repetition subsample is drawn in the builder below.
  ## "levels": keep exactly the observations at the levels in S. A retained
  ## cell observes every level in S (the mask) and each at most once (the
  ## check above), so it keeps exactly L0 observations, and level s is the
  ## same level in every cell.
  idx <- if (trim == "random") which(in_block[cell] > 0L)
         else which(in_block[cell] > 0L & lev %in% S)
  if (!length(idx))
    stop("Internal error: the selected blocks contain no observations.",
         call. = FALSE)   # nocov
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
  ## cell clears it, so this is the identity, but the group order below must
  ## never be set by a block with nothing in it).
  used_q <- sort(unique(blk_k))
  blocks <- blocks[used_q]
  blk_k <- match(blk_k, used_q)
  min_side <- min(vapply(blocks,
                         function(b) min(length(b$rows), length(b$cols)),
                         integer(1)))

  if (trim == "random") {
    cell_d <- .dense_id(cell_k)                 # dense id among retained cells
    ncell_d <- max(cell_d)
    ell_d <- tabulate(cell_d, nbins = ncell_d)  # all >= L0 by the mask
    key_k <- ord_key[idx]
    ## The block builder's stride, so that sub-seed offset 1 here can never
    ## coincide with a block seed (offsets 4q - 1, 4q) of another repetition.
    stride <- max(1000, 4 * length(blocks) + 1)
    perm_builder <- function(rep_seed, K) {
      take <- .downsample_to_L0(cell_d, ell_d, L0,
                                seed = .sub_seed(rep_seed, 1L, stride))
      ## the slot: rank of each survivor inside its cell, by `rep` order
      slot <- .within_cell_slot(cell_d[take], key_k[take], ncell_d)
      ops <- .build_obs_perms_blocks(rep_seed, K, blocks,
                                     ri = ri_k[take], ci = ci_k[take],
                                     blk = blk_k[take],
                                     lrow = lrow_k[take], lcol = lcol_k[take],
                                     permute = "both", slot = slot)
      attr(ops, "rows") <- take             # this repetition's rows of idx
      ops
    }
  } else {
    slot <- match(lev[idx], S)                 # 1..L0, in level order
    perm_builder <- function(rep_seed, K)
      .build_obs_perms_blocks(rep_seed, K, blocks,
                              ri = ri_k, ci = ci_k, blk = blk_k,
                              lrow = lrow_k, lcol = lcol_k, permute = "both",
                              slot = slot)
  }

  list(idx = idx, S = S,
       S_labels = if (trim == "levels") S_labels_all[S] else NULL,
       cell = cell, ncell = ncell, ell = ell, n_mask = sum(mask),
       in_block = in_block, blocks = blocks, blk = blk_k,
       min_side = min_side, n_cells_used = n_cells_used,
       n_row = max(ri), n_col = max(ci), trim = trim,
       perm_builder = perm_builder)
}
