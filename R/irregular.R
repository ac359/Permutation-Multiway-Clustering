#' Invariant permutation test for irregular two-way layouts (Section 6.4)
#'
#' Finite-sample valid test of \eqn{H_0: \beta = b} for a two-way layout
#' \deqn{y_{ijl} = x_{ijl}^\top \gamma + d_{ijl}^\top \beta + \varepsilon_{ijl}}
#' where cell (i, j) holds ell_ij observations indexed by l, the cell sizes are
#' unequal, and permuting \emph{within} a cell is either invalid or powerless.
#' This is the procedure of Guo, Toulis and Wang (2026), Section 6.4. Two cases
#' motivate it, and \code{\link{mwperm_layout}} handles neither:
#' \itemize{
#'   \item the replication index l is really \emph{time}, so the errors are not
#'     exchangeable across l (invariance InvB, not InvA) and within-cell
#'     permutation is \strong{invalid};
#'   \item d is constant within each (i, j) cell (a dyad-level covariate), so
#'     within-cell permutation leaves the residual statistic unchanged and the
#'     test has \strong{no power}.
#' }
#' In both cases the fix is to permute the \emph{cells}, across i and j,
#' exactly as \code{\link{mwperm_dyadic}} does -- which needs equal cell sizes
#' and a complete array, neither of which an irregular layout has. Section 6.4
#' obtains both by combining the missing-data machinery of Procedure 2 with the
#' panel construction:
#' \enumerate{
#'   \item form the cell sizes ell_ij and, for a threshold L0, the mask
#'     M_ij = 1 if ell_ij >= L0, and 0 otherwise;
#'   \item run the biclique search (\code{\link{find_bicliques}}, Algorithm 2)
#'     on that mask to obtain disjoint fully observed blocks I_q by J_q;
#'   \item inside each selected block, discard ell_ij - L0 observations from
#'     each cell uniformly at random, leaving exactly L0 everywhere;
#'   \item apply Procedure 2 to what remains: a row group on I_q and a column
#'     group on J_q per block, of common order K + 1, applied identically across
#'     the L0 within-cell slots, so cell (i, j) slot l maps to cell
#'     (pi(i), sigma(j)) slot l.
#' }
#' Holding the slot fixed is the same device \code{\link{mwperm_panel}} uses
#' for time, and it is what makes the test valid when l indexes periods.
#'
#' @section Assumptions:
#' This test does \strong{not} assume within-cell exchangeability. What it needs
#' is exchangeability across the cell indices (i, j) \emph{within each slot}:
#' for every slot l, the errors in a block satisfy
#'
#'   (eps_ijl for i in I_q, j in J_q)  has the same distribution, given X and D,
#'   as (eps_[pi(i)][sigma(j)]l for i in I_q, j in J_q)
#'
#' with the \emph{same} (pi, sigma) used in every slot -- invariance (InvB) of
#' Section 6.2, applied blockwise. That is required together with Assumption 4
#' on the mask: the missingness pattern (here, which cells clear the L0
#' threshold) is independent of the errors given the covariates. Under those two
#' conditions the p-value is exact in finite samples by the argument of Theorem
#' 4. Arbitrary slot effects shared across cells -- an arbitrary common time
#' trend -- are permitted, which is the whole point.
#'
#' @section Choosing L0:
#' L0 trades cells against within-cell depth: a small L0 keeps more cells in the
#' mask (so larger blocks, larger K, finer p-value resolution) but throws away
#' more observations per cell; a large L0 keeps deeper cells but fewer of
#' them. The paper recommends tuning it by grid
#' search to minimise the loss of observations. There is no free lunch in
#' choosing it from the data on the \emph{outcome}; select it from the cell
#' sizes alone, which are ancillary under Assumption 4. The \code{note} field of
#' the fitted object reports how many cells and observations survived.
#'
#' @section Reproducibility:
#' Step (iii) deletes observations at random. That draw is taken \strong{once},
#' from \code{seed}, before the permutation loop, so the \code{n_reps}
#' repetitions average over the permutation-group draw but not over the
#' subsample. To average over the subsample as well -- as Section B of the paper
#' does -- run the test under several seeds and take the median p-value across
#' those runs.
#'
#' @inheritParams mwperm_dyadic
#' @param d Numeric vector or matrix of the covariate(s) of interest. Unlike
#'   \code{\link{mwperm_layout}}, \code{d} \emph{may} be constant within cells;
#'   that is one of the cases this design exists for, and no warning is issued.
#' @param x Optional numeric matrix or data frame of nuisance covariates; an
#'   intercept is always added internally. May be \code{NULL}.
#' @param row,col Cell identifiers along the two layout dimensions.
#' @param rep Optional within-cell index (the replication or period
#'   identifier). It fixes which slot each observation occupies, and hence
#'   which observations are aligned across cells; when \code{NULL}, order of
#'   appearance within the cell is used. Supply it when the within-cell index
#'   means something -- a period, a wave -- so that slot l is the same period in
#'   every cell.
#' @param L0 Integer threshold, at least 2: cells with fewer than \code{L0}
#'   observations are masked out, and every retained cell is reduced to exactly
#'   \code{L0}. Required.
#' @param K Number of non-identity permutations; defaults to the smallest block
#'   side over the selected blocks -- that is, the smallest of |I_q| and |J_q|
#'   over all q -- minus one, capped at 199. Must satisfy \code{K + 1 <=} that
#'   smallest block side.
#' @param min_block Minimum block side(s) for the biclique search; see
#'   \code{\link{find_bicliques}}.
#' @param block_method \code{"greedy"} (default) or \code{"exact"}; see
#'   \code{\link{find_bicliques}}.
#'
#' @return An object of class \code{"mwperm"}, with the extra fields
#'   \code{n_blocks}, \code{cells_used}, \code{cells_total} and \code{L0}.
#'   \code{estimate}/\code{se_naive} are the OLS estimate and naive SE on the
#'   retained data, \code{conf_int}/\code{conf_set} (or
#'   \code{conf_region}/\code{conf_box} for several coefficients) the IPT
#'   inverted-test confidence set, and \code{pvalue} the IPT permutation
#'   p-value; see \code{\link{mwperm_dyadic}} for the field provenance in full.
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation
#'   inference under multi-way clustering and missing data, Section 6.4 and
#'   Procedure 2. arXiv:2601.08610.
#' @seealso \code{\link{mwperm_layout}} (Section 6.3, within-cell permutation),
#'   \code{\link{mwperm_missing}} (Procedure 2 with one observation per cell),
#'   \code{\link{mwperm_panel}} (Section 6.2, complete balanced panels).
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
    stop(paste0("`L0` is required: it is the cell-size threshold that ",
                "defines the mask M_ij = 1{ell_ij >= L0} (Section 6.4). ",
                "Pick it from the cell sizes -- see the `Choosing L0` ",
                "section of ?mwperm_irregular."), call. = FALSE)
  L0 <- suppressWarnings(as.integer(L0))
  if (length(L0) != 1L || is.na(L0) || L0 < 2L)
    stop("`L0` must be a single integer >= 2.", call. = FALSE)

  ri <- .dense_id(row, "row")
  ci <- .dense_id(col, "col")
  n_row <- max(ri)
  n_col <- max(ci)

  ## --- (i) cell sizes and the mask M_ij = 1{ell_ij >= L0} --------------------
  cell <- .dense_id(interaction(ri, ci, drop = TRUE))   # dense occupied-cell id
  ncell <- max(cell)
  ell <- tabulate(cell, nbins = ncell)                  # cell sizes
  first <- match(seq_len(ncell), cell)                  # a row of each cell
  mask <- ell >= L0                                     # M_ij
  if (!any(mask))
    stop(sprintf(paste0("No cell has at least L0 = %d observations (the ",
                        "largest cell has %d), so the mask M_ij = 1{ell_ij ",
                        ">= L0} is empty. Lower `L0`."),
                 L0, max(ell)), call. = FALSE)

  ## --- (ii) biclique search on the mask (Algorithm 2) ------------------------
  mrow <- ri[first[mask]]
  mcol <- ci[first[mask]]                               # cells that clear L0
  blocks <- find_bicliques(mrow, mcol, min_block = min_block,
                           method = block_method)
  if (length(blocks) == 0L)
    stop(sprintf(paste0("No fully observed block with both sides >= ",
                        "min_block = %s was found among the %d cells that ",
                        "clear L0 = %d. Lower `min_block`, lower `L0` (more ",
                        "cells clear the threshold), or use ",
                        "mwperm_layout() if within-cell permutation is valid ",
                        "for these data."),
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
  keep_cell <- in_block > 0L
  sel <- which(keep_cell[cell])              # observations in a selected cell
  if (!length(sel))
    stop("Internal error: the selected blocks contain no observations.",
         call. = FALSE)   # nocov

  ## --- (iii) reduce every retained cell to exactly L0 observations -----------
  ## Uniform, reproducible from `seed`, and with the same RNG-state hygiene as
  ## .downsample_to_L0() (which is what does the drawing): the caller's random
  ## stream is left exactly as it was found.
  cell_sel <- .dense_id(cell[sel])
  ell_sel <- tabulate(cell_sel, nbins = max(cell_sel))
  take <- .downsample_to_L0(cell_sel, ell_sel, L0, seed = seed)
  idx <- sel[take]                           # retained observations, sorted

  yk <- y[idx]
  Dk <- D[idx, , drop = FALSE]
  Xk <- X[idx, , drop = FALSE]
  ri_k <- ri[idx]
  ci_k <- ci[idx]
  cell_k <- .dense_id(cell[idx])
  Nk <- length(idx)
  n_cells_used <- max(cell_k)

  ## --- (iv) slot index, held fixed by the permutation ------------------------
  ## Slot l of cell (i, j) maps to slot l of cell (pi(i), sigma(j)); `rep`
  ## decides which observation is which slot, so that slot l means the same
  ## period in every cell when `rep` is a period identifier.
  slot_k <- .within_cell_slot(cell_k, if (is.null(rep)) NULL else rep[idx],
                              n_cells_used)

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

  ## --- group order: the smallest permuted block side over the blocks ---------
  used_q <- sort(unique(blk_k))
  blocks <- blocks[used_q]                   # blocks that actually kept cells
  blk_k <- match(blk_k, used_q)
  min_side <- min(vapply(blocks,
                         function(b) min(length(b$rows), length(b$cols)),
                         integer(1)))
  K_was_null <- is.null(K)
  K <- .default_K(K, min_side)

  sizes <- vapply(blocks, function(b) c(length(b$rows), length(b$cols)),
                  integer(2))
  note <- c(
    sprintf(paste0("Section 6.4 (irregular design), L0 = %d: %d of %d cells ",
                   "clear the threshold, %d of those lie in the %d fully ",
                   "observed block%s the biclique search extracted, and each ",
                   "was reduced to exactly %d observations. Kept %d of %d ",
                   "observations (%.1f%%). The discarded observations are ",
                   "what buys exact validity under an unequal design."),
            L0, sum(mask), ncell, n_cells_used, length(blocks),
            if (length(blocks) == 1L) "" else "s", L0, Nk, N, 100 * Nk / N),
    sprintf("Block sizes (rows x cols): %s.",
            paste(sprintf("%dx%d", sizes[1, ], sizes[2, ]), collapse = ", ")),
    if (K_was_null && is.numeric(alpha) && length(alpha) == 1L &&
        is.finite(alpha) && alpha > 0 && alpha < 1 && 1 / (K + 1) > alpha)
      sprintf(paste0(
        "The smallest selected block (permuted side %d) caps the group ",
        "order: K = %d, so no rejection is attainable at alpha = %.3g. A ",
        "larger L0 gives deeper cells but fewer of them; a smaller L0 admits ",
        "more cells and can support larger blocks. Raise `min_block` to stop ",
        "small blocks from setting K."),
        min_side, K, alpha)
    else character(0))

  ## --- (v) Procedure 2 on the stacked retained data --------------------------
  perm_builder <- function(rep_seed)
    .build_obs_perms_blocks(rep_seed, K, blocks,
                            ri = ri_k, ci = ci_k, blk = blk_k,
                            lrow = lrow_k, lcol = lcol_k, permute = "both",
                            slot = slot_k)

  res <- .ipt_engine(yk, Dk, Xk, perm_builder, K = K, n_reps = n_reps,
                     seed = seed, alpha = alpha, conf_int = conf_int,
                     beta_null = beta_null, grid = grid,
                     type = "irregular (Section 6.4)", d_names = d_names,
                     n_clusters = c(row = n_row, col = n_col), call = cl,
                     n_cores = n_cores, ci_agg = aggregate)
  res$note <- c(note, res$note)
  res$n_blocks <- length(blocks)
  res$cells_used <- n_cells_used
  res$cells_total <- ncell
  res$L0 <- L0
  res
}
