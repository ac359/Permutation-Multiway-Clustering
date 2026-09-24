## ============================================================================
## R/irregular.R -- irregular layouts: Section 6.4, step (i) as printed
##
## Purpose. mwperm_irregular() handles unequal cell sizes when permuting
##   within cells is powerless (d constant within a cell): .irregular_design()
##   masks the cells with at least L0 observations, finds fully observed
##   blocks, and returns a builder that, in EVERY repetition, cuts each
##   retained cell to L0 observations at random and builds the block-diagonal
##   group over the survivors with their within-cell rank held fixed. The
##   builder names the survivors in attr(, "rows"), which is how the engine
##   runs that repetition on the subsample.
## Paper. Guo, Toulis & Wang (2026), Section 6.4: the mask
##   M_ij = 1{ell_ij >= L0}, step (i) (Algorithm 2 and the random drop of
##   ell_ij - L0 observations), step (ii) (Procedure 2 on what remains), and
##   the median over repeated runs (Remark 1; 100 runs in their Appendix B).
##   Valid by the argument of Theorem 4 when the observations inside a cell
##   are exchangeable replicates; periods with a common effect belong to
##   mwperm_panel_missing(time =, L0 =) instead (each fit's first note says so).
## Pipeline. mwperm() -> dispatch -> [design worker] -> [permutation
##   construction, via .build_obs_perms_blocks() in missing.R] -> projection
##   engine -> median aggregation -> test inversion -> S3 methods.
## ============================================================================

#' Invariant permutation test for irregular two-way layouts (Section 6.4)
#'
#' Finite-sample valid test of H0: beta = b for a two-way layout
#' \preformatted{  y_ijl = x_ijl' gamma + d_ijl' beta + eps_ijl}
#'
#' where cell (i, j) holds ell_ij observations indexed by l, the cell sizes
#' are unequal, and permuting *within* a cell is powerless -- typically
#' because d is constant within each (i, j) cell (a dyad-level covariate),
#' so within-cell permutation leaves the residual statistic unchanged and
#' [mwperm_layout()] has **no power**. This is the procedure of Guo, Toulis
#' and Wang (2026), Section 6.4, step (i) exactly as printed and as the paper
#' applies it in its Appendix B.
#'
#' The fix is to permute the *cells*, across i and j, exactly as
#' [mwperm_dyadic()] does -- which needs equal cell sizes and a complete
#' array, neither of which an irregular layout has. Section 6.4 obtains both
#' by combining the missing-data machinery of Procedure 2 with the panel
#' construction of Case (B):
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
#' *Aggregation over repetitions* in [mwperm_dyadic()]). The paper's Appendix
#' B used 100 repetitions; the package default is `n_reps = 500`, and a
#' final run is worth 1000 (runtime is linear in `n_reps`).
#'
#' @details Implements Section 6.4 of Guo, Toulis and Wang (2026). Step (i):
#'   the mask `M_ij = 1{ell_ij >= L0}`, Algorithm 2, and a random cut of every
#'   retained cell to L0 observations, redrawn in every repetition. Step (ii):
#'   Procedure 2 with the within-cell position held fixed. The median over
#'   repetitions follows their Remark 1.
#' @section Assumptions: This test does **not** assume within-cell
#'   exchangeability. What it needs is exchangeability of the retained error
#'   array across the cell indices (i, j) within each block, position by
#'   position: with the *same* (pi, sigma) at every within-cell position l,
#'
#' (eps_ijl for i in I_q, j in J_q) has the same distribution, given X and D,
#' as (eps_[pi(i)][sigma(j)]l for i in I_q, j in J_q),
#'
#' together with Assumption 4 on the mask: which cells clear L0 is
#' independent of the errors given the covariates. Under those two conditions
#' the p-value is exact in finite samples by the argument of Theorem 4. The
#' position l of a retained observation is its rank among the survivors of a
#' random draw, so the condition holds when the observations inside a cell
#' are **exchangeable replicates** -- individuals sampled within a cell, as
#' in the paper's Appendix B, with no effect attached to the within-cell
#' index -- with any cell-level dependence of the random-effects kind
#' (eps_ijl = eta_i + xi_j + u_ijl).
#'
#' @section When the repeats are periods: When the within-cell index is a
#'   time period, or anything else with an effect zeta_l shared across
#'   cells, the condition fails: the random trim keeps different periods in
#'   different cells, the position it holds fixed is no longer the period,
#'   and the test can over-reject **even with a cell-constant `d`**. In this
#'   package's check (30 x 30 cells, rows 1-15 observed in periods 1-5 and
#'   rows 16-30 in periods 4-8, a period effect of 0.5 per period, a
#'   cell-constant `d` correlated with the row cohort, `L0 = 3`) the size at
#'   nominal 0.05 was 0.39, 0.375 at the default `n_reps = 500`, and 0.105
#'   with a period effect of only 0.2; [mwperm_panel_missing()] on the same
#'   data gave 0.028. A `d` that also varies over the periods (staggered
#'   adoption) is worse still: size 1.000, against 0.040 for the aligned cut.
#'   Earlier checks that found the trim at nominal under a strong period
#'   effect (0.043 with a cell-constant `d`, 0.051 with an i.i.d. one) held
#'   only because every cell observed the same periods. The fit cannot see a
#'   period effect, so its note states this condition on every fit. For any
#'   period index use [mwperm_panel_missing()] with the period as `time` and
#'   `L0 =` the number of periods to keep: it retains the *same* L0 periods
#'   in every cell (chosen from the observation pattern alone) and holds the
#'   period fixed, which is condition InvB. (Before 0.4.2 that cut was this
#'   function's `trim = "levels"` option; `mwperm_panel_missing(L0 = ,
#'   time_fe = FALSE)` reproduces it exactly.)
#'
#' @section Choosing L0: L0 trades cells against within-cell depth: a small L0
#'   keeps more cells in the mask (so larger blocks, larger K, finer p-value
#'   resolution) but throws away more observations per cell; a large L0 keeps
#'   deeper cells but fewer of them. The paper recommends tuning it by grid
#'   search to minimise the loss of observations. There is no free lunch in
#'   choosing it from the data on the *outcome*; select it from the cell
#'   sizes alone, which are ancillary under Assumption 4. The `note` field of
#'   the fitted object reports how many cells and observations survived.
#'
#' @section Reproducibility and the confidence set: The subsample of each
#'   repetition is drawn from that repetition's seed (rep r uses
#'   `seed + r - 1`), with the usual RNG-state hygiene: a seeded call leaves
#'   the caller's random stream exactly as it found it, and the rep-parallel
#'   path (`n_cores > 1` with a `seed`) reproduces the serial result bit for
#'   bit. At the default `n_reps = 500` the confidence interval is found by
#'   outward bracketing and bisection unless `K` is small (`K <= 14`): the
#'   exact set evaluates about `2 K^2 n_reps` candidate end points, which
#'   exceeds the engine's budget beyond that, and its cost grows like
#'   `K^3 n_reps^2`. The fit records which path it took in `ci_method` and
#'   in a note; the end points agree to the bisection tolerance.
#'
#' @inheritParams mwperm_dyadic
#' @param d Numeric vector or matrix of the covariate(s) of interest. Unlike
#'   [mwperm_layout()], `d` *may* be constant within cells; that is one of the
#'   cases this design exists for, and no warning is issued.
#' @param x Optional numeric matrix or data frame of nuisance covariates; an
#'   intercept is always added internally. May be `NULL`.
#' @param row,col Cell identifiers along the two layout dimensions.
#' @param rep Optional within-cell index (the replication identifier). It
#'   only orders the observations inside a cell (the survivors of each
#'   repetition's draw are ranked by it; order of appearance when `NULL`), so
#'   the fit does not depend on it beyond that order.
#' @param L0 Integer, at least 2: the number of observations retained per
#'   cell. Cells with fewer than `L0` observations are masked out, and every
#'   retained cell is reduced to exactly `L0`. Required.
#' @param K Number of non-identity permutations; defaults to the smallest
#'   block side over the selected blocks -- that is, the smallest of |I_q| and
#'   |J_q| over all q -- minus one, capped at 199. Must satisfy `K + 1 <=`
#'   that smallest block side.
#' @param min_block Minimum block side(s) for the biclique search; see
#'   [find_bicliques()].
#' @param block_method `"greedy"` (default) or `"exact"`; see
#'   [find_bicliques()].
#' @param n_reps Number of repetitions, each with its own random per-cell
#'   trim and its own permutation group; the reported p-value is the median
#'   over them (Remark 1 of the paper) and the confidence set inverts that
#'   median. Defaults to 500 -- the paper's Appendix B used 100, which the
#'   package's advisors consider the low end; use 1000 for a final run. Every
#'   other front end defaults to 10, because only here is the data itself
#'   redrawn in every repetition.
#'
#' @return An object of class `"mwperm"`, with the extra fields `n_blocks`,
#'   `cells_used`, `cells_total` and `L0`. `n_obs` is the number of
#'   observations each repetition is computed on, `cells_used * L0`.
#'   `estimate`/`se_naive` are the OLS estimate and naive SE on *all*
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
#'   [mwperm_panel_missing()] (the same L0 periods in every cell, for
#'   repeats that are periods), [mwperm_panel()] (Section 6.2, complete
#'   balanced panels).
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
                             n_reps = 500L, seed = NULL, grid = NULL,
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
    stop(paste0("`L0` is required: it is the cell-size threshold that ",
                "defines the mask M_ij = 1{ell_ij >= L0} and the number of ",
                "observations every retained cell is reduced to (Section ",
                "6.4). Pick it from the cell sizes -- see the `Choosing L0` ",
                "section of ?mwperm_irregular."), call. = FALSE)
  L0 <- suppressWarnings(as.integer(L0))
  if (length(L0) != 1L || is.na(L0) || L0 < 2L)
    stop("`L0` must be a single integer >= 2.", call. = FALSE)

  ## --- (i)-(iii) the mask, the blocks and the cut ----------------------------
  pd <- .irregular_design(row, col, rep, L0, min_block = min_block,
                          block_method = block_method)
  idx <- pd$idx                              # observations of the retained cells
  yk <- y[idx]
  Dk <- D[idx, , drop = FALSE]
  Xk <- X[idx, , drop = FALSE]
  blocks <- pd$blocks
  n_cells_used <- pd$n_cells_used
  Nk <- n_cells_used * L0                    # observations per repetition
  ## The engine checks N > 2p on the data it is handed, which is every
  ## observation of the retained cells; the projection in each repetition
  ## sees only cells_used x L0 of them, so that is the count that must clear
  ## 2p.
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

  ## Does `d` vary inside any retained cell? The same diagnostic
  ## mwperm_layout() computes for its no-power warning, read the other way
  ## round: there a cell-constant d is the problem, here it is the case the
  ## design exists for. The random trim is exact when the observations inside
  ## a cell are exchangeable replicates; if the within-cell index is a period
  ## with an effect shared across cells the test can over-reject even with a
  ## cell-constant d (the first note says so on every fit), and a d that
  ## varies over the periods makes that far worse. Only d varying is visible
  ## in the data, so it is reported as a note, never a warning.
  cell_k <- pd$cell[idx]
  within_var <- tapply(seq_along(idx), cell_k, function(ii)
    any(apply(Dk[ii, , drop = FALSE], 2L, function(v) diff(range(v)) > 0)))
  d_varies <- any(unlist(within_var))

  sizes <- vapply(blocks, function(b) c(length(b$rows), length(b$cols)),
                  integer(2))
  note <- c(
    sprintf(paste0("Section 6.4 (irregular design), L0 = %d: %d of %d ",
                   "cells have at least L0 observations, %d of those lie ",
                   "in the %d fully observed block%s the biclique search ",
                   "extracted, and in every repetition each of them is cut ",
                   "to exactly %d observations drawn at random (a fresh ",
                   "draw per repetition; the reported p-value is the ",
                   "median over %d). Each repetition uses %d of the %d ",
                   "observations (%.1f%%); the OLS estimate uses all %d ",
                   "observations of the retained cells. The random trim is ",
                   "exact only when the observations inside a cell are ",
                   "exchangeable replicates: if the within-cell index is a ",
                   "period, or anything else with an effect shared across ",
                   "cells, the test can over-reject even with a ",
                   "cell-constant `d` -- use mwperm_panel_missing(time = ",
                   "<rep>, L0 = <L0>), which keeps the same L0 periods in ",
                   "every cell."),
            L0, pd$n_mask, pd$ncell, n_cells_used, length(blocks),
            if (length(blocks) == 1L) "" else "s", L0, n_reps, Nk, N,
            100 * Nk / N, length(idx)),
    sprintf("Block sizes (rows x cols): %s.",
            paste(sprintf("%dx%d", sizes[1, ], sizes[2, ]), collapse = ", ")),
    if (d_varies)
      paste0("`d` varies within cells. If the within-cell index is a ",
             "period, or anything else with an effect shared across cells, ",
             "this test is NOT valid, and a `d` that varies over the ",
             "periods makes its over-rejection severe; use ",
             "mwperm_panel_missing(time = <rep>, L0 = <L0>), which keeps ",
             "the same L0 periods in every cell.")
    else character(0),
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
  res$n_obs <- Nk                            # per repetition
  res$n_blocks <- length(blocks)
  res$cells_used <- n_cells_used
  res$cells_total <- pd$ncell
  res$L0 <- L0
  res
}

#' Steps (i)-(iii) of Section 6.4: the mask, the blocks, the cut.
#'
#' Everything `mwperm_irregular()` does to the observation pattern before the
#' engine sees a number, factored out so that the structural claims -- every
#' repetition keeps exactly L0 observations per retained cell, and every
#' gather vector preserves the within-cell position -- can be asserted on the
#' front end's OWN retained set and permutation builder rather than on a
#' hand-built imitation of them (see `tests/test-irregular.R`).
#'
#' The mask is the printed count mask 1{ell_ij >= L0}, `idx` is EVERY
#' observation of the retained cells, and the cut happens inside
#' `perm_builder(rep_seed, K)`: it draws that repetition's subsample
#' (`.downsample_to_L0()`, seeded from the rep seed at sub-seed offset 1 --
#' the block builder uses offsets 4q - 1 and 4q, so nothing collides), ranks
#' the survivors inside each cell by `rep` (order of appearance when `NULL`),
#' builds the block-diagonal Procedure 2 group over the survivors with that
#' rank as the slot held fixed, and returns it with the attribute `"rows"`
#' naming the survivors within `idx`, which is how the engine learns to
#' compute that repetition on the subsample. (The deterministic level-aligned
#' cut this function carried as `trim = "levels"` until 0.4.1 is now
#' `.panel_missing_design(L0 = )` in `R/panel_missing.R`.)
#'
#' @details Section 6.4, step (i) of GTW (2026), and the Procedure 2 group of
#'   step (ii).
#' @param row,col cell identifiers, one per observation.
#' @param rep within-cell order identifier, or `NULL` for order of appearance.
#' @param L0 number of observations retained per cell.
#' @param min_block,block_method passed to `find_bicliques()`.
#' @return a list: `idx` (sorted indices of the observations handed to the
#'   engine), `cell`, `ncell`, `ell` (observations per cell), `n_mask`,
#'   `in_block` (block of each cell, 0 = discarded), `blocks` (those that
#'   retained cells), `min_side`, `n_cells_used`, `n_row`, `n_col`, and
#'   `perm_builder(rep_seed, K)`, the Procedure 2 builder described above.
#' @keywords internal
#' @noRd
.irregular_design <- function(row, col, rep = NULL, L0, min_block = 3L,
                              block_method = "greedy") {
  ri <- .dense_id(row, "row")
  ci <- .dense_id(col, "col")
  N <- length(ri)
  cell <- .dense_id(interaction(ri, ci, drop = TRUE))   # dense occupied-cell id
  ncell <- max(cell)

  ## --- (i) the mask ----------------------------------------------------------
  ## The printed mask: cell (i, j) clears the threshold when it holds at
  ## least L0 observations, whatever their labels. `rep` only orders the
  ## observations inside a cell (the survivors are ranked by it).
  ell <- tabulate(cell, nbins = ncell)                  # observations per cell
  if (max(ell) < L0)
    stop(sprintf(paste0("No cell has at least L0 = %d observations (the ",
                        "largest cell has %d), so the mask M_ij = 1{ell_ij ",
                        ">= L0} is empty. Lower `L0`."),
                 L0, max(ell)), call. = FALSE)
  mask <- ell >= L0                  # Section 6.4: M_ij = 1{ell_ij >= L0}
  ord_key <- if (is.null(rep)) seq_len(N) else as.numeric(factor(rep))
  first <- match(seq_len(ncell), cell)                  # a row of each cell

  ## --- (ii) biclique search on the mask (Algorithm 2) ------------------------
  blocks <- find_bicliques(ri[first[mask]], ci[first[mask]],
                           min_block = min_block, method = block_method)
  if (length(blocks) == 0L)
    stop(sprintf(paste0("No fully observed block with both sides >= ",
                        "min_block = %s was found among the %d cells that ",
                        "clear L0 = %d. Lower `min_block`, lower `L0` (more ",
                        "cells clear the threshold), or use mwperm_layout() ",
                        "if within-cell permutation is valid for these ",
                        "data."),
                 paste(min_block, collapse = "x"), sum(mask), L0),
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
  ## Every observation of a retained cell is handed to the engine and the
  ## per-repetition subsample is drawn in the builder below.
  idx <- which(in_block[cell] > 0L)
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

  cell_d <- .dense_id(cell_k)                 # dense id among retained cells
  ncell_d <- max(cell_d)
  ell_d <- tabulate(cell_d, nbins = ncell_d)  # all >= L0 by the mask
  key_k <- ord_key[idx]
  ## The block builder's stride, so that sub-seed offset 1 here can never
  ## coincide with a block seed (offsets 4q - 1, 4q) of another repetition.
  stride <- max(1000, 4 * length(blocks) + 1)
  ## Step (i)'s random cut, redrawn per repetition from that repetition's
  ## seed, then step (ii)'s Procedure 2 group over the survivors, with the
  ## within-cell rank as the slot held fixed: survivor l of cell (i, j) maps
  ## to survivor l of cell (pi(i), sigma(j)).
  perm_builder <- function(rep_seed, K) {
    take <- .downsample_to_L0(cell_d, ell_d, L0,
                              seed = .sub_seed(rep_seed, 1L, stride))
    ## the slot: rank of each survivor inside its cell, by `rep` order
    slot <- .within_cell_slot(cell_d[take], key_k[take], ncell_d,
                              ranked = TRUE)   # key_k is already a rank
    ops <- .build_obs_perms_blocks(rep_seed, K, blocks,
                                   ri = ri_k[take], ci = ci_k[take],
                                   blk = blk_k[take],
                                   lrow = lrow_k[take], lcol = lcol_k[take],
                                   permute = "both", slot = slot)
    attr(ops, "rows") <- take             # this repetition's rows of idx
    ops
  }

  list(idx = idx, cell = cell, ncell = ncell, ell = ell, n_mask = sum(mask),
       in_block = in_block, blocks = blocks, blk = blk_k,
       min_side = min_side, n_cells_used = n_cells_used,
       n_row = max(ri), n_col = max(ci),
       perm_builder = perm_builder)
}
