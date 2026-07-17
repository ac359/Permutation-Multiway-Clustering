#' Invariant permutation test for dyadic regression with missing cells
#'
#' Finite-sample valid test of \eqn{H_0: \beta = b} in the dyadic model when the
#' design is \emph{not} a complete array, i.e. some \eqn{(i, j)} cells are
#' unobserved (Procedure 2 of Guo, Toulis and Wang, 2026). General row/column
#' permutations are infeasible once cells are missing, so the test is built on
#' \emph{fully observed bicliques}: maximal-ish blocks of row clusters
#' \eqn{I_q} and column clusters \eqn{J_q} for which every cell
#' \eqn{I_q \times J_q} is observed. The block-cyclic permutation group
#' (\code{\link{build_perm_set}}) is applied within each block, the blocks are
#' chosen to be disjoint in both rows and columns so that the joint permutation
#' is a genuine relabelling of exchangeable clusters, and the residual
#' statistics are pooled across blocks. Cells outside every selected block are
#' discarded; this is the price of validity under missingness and is reported in
#' the result.
#'
#' Finding the largest fully observed biclique is NP-hard, so a greedy
#' heuristic is used; an approximate (sub-maximal) solution still yields a valid
#' test, only a less powerful one. \strong{The smallest selected block is the
#' binding constraint on resolution}: the group order \code{K + 1} is capped at
#' its smaller side, so the smallest attainable p-value is
#' \code{1/(K + 1) >= 1/min_side} no matter how large the other blocks are.
#' Adding small blocks pools more data yet can destroy resolution -- rejecting
#' at level \code{alpha} requires a fully observed block with both sides at
#' least \code{ceiling(1/alpha)} (20 for \code{alpha = 0.05}). When the
#' default \code{K} makes rejection at \code{alpha} unattainable the fit says
#' so in its \code{note}; raise \code{min_block} to stop small blocks from
#' setting \code{K}, at the price of discarding their cells.
#'
#' @inheritParams mwperm_dyadic
#' @param row,col Cluster identities of each observed cell. Cells that are absent
#'   from the data are treated as missing.
#' @param min_block Integer; only bicliques whose smaller side is at least
#'   \code{min_block} are used. Larger values give finer p-value resolution
#'   (larger \code{K}) but discard more data. Defaults to 3.
#' @param block_method Passed to \code{\link{find_bicliques}}: \code{"greedy"}
#'   (default) for the fast heuristic, or \code{"exact"} for a branch-and-bound
#'   search of maximum-area blocks (with a node-budget fallback to greedy).
#' @param K Number of non-identity permutations. Defaults to
#'   \code{min_q min(|I_q|, |J_q|) - 1} over the selected blocks, capped at 199.
#'   Must not exceed that quantity.
#'
#' @return An object of class \code{"mwperm"}: \code{estimate}/\code{se_naive}
#'   are the OLS estimate and naive SE, \code{conf_int} (or
#'   \code{conf_region}/\code{conf_box} for several coefficients) the IPT
#'   inverted-test confidence set, and \code{pvalue} the IPT permutation
#'   p-value; see \code{\link{mwperm_dyadic}} for the field provenance in
#'   full. Its \code{n_obs} is the number of cells actually used (inside the
#'   selected blocks), and the \code{note} field records how many cells and
#'   blocks were retained.
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation
#'   inference under multi-way clustering and missing data. arXiv:2601.08610.
#'
#' @seealso \code{\link{mwperm_dyadic}} for the complete-array case;
#'   \code{\link{find_bicliques}} for the block finder used internally.
#'
#' @examples
#' data(trade_dyadic)
#' ## drop self-trade and a few random links to create missingness
#' set.seed(1)
#' d <- trade_dyadic[trade_dyadic$importer != trade_dyadic$exporter, ]
#' d <- d[sample(nrow(d), round(0.9 * nrow(d))), ]
#' fit <- with(d, mwperm_missing(y = log_trade, d = log_dist,
#'                               x = cbind(log_gdp_i, log_gdp_j),
#'                               row = importer, col = exporter,
#'                               min_block = 3, seed = 1))
#' fit
#' @export
mwperm_missing <- function(y, d, x = NULL, row, col, K = NULL,
                           alpha = 0.05, beta_null = 0, conf_int = TRUE,
                           n_reps = 10L, seed = NULL, grid = NULL,
                           min_block = 3L, block_method = c("greedy", "exact"),
                           n_cores = 1L) {
  cl <- match.call()
  block_method <- match.arg(block_method)
  y <- .check_y(y); N <- length(y)
  D <- as.matrix(d); d_names <- .coef_names(D, deparse(substitute(d)))
  X <- .make_X(x, N)
  .check_lengths(N, list(row = row, col = col))
  ## Validate the FULL supplied data before the biclique step discards cells
  ## (audit F5.4): the observed cells are the user's data contract, and an NA
  ## in a soon-to-be-discarded cell should not pass silently.
  .check_finite(list(y = y, d = D, x = X))

  ri <- .dense_id(row, "row"); ci <- .dense_id(col, "col")
  if (anyDuplicated(cbind(ri, ci)))
    stop("Expected at most one observation per (row, col) cell.", call. = FALSE)
  n_row <- max(ri); n_col <- max(ci)

  ## --- find disjoint fully observed bicliques --------------------------------
  blocks <- find_bicliques(ri, ci, min_block = min_block, method = block_method)
  if (length(blocks) == 0L)
    stop(sprintf(paste0("No fully observed block with both sides >= min_block = %d ",
                        "was found. Lower `min_block` (the test needs at least a ",
                        "2x2 observed block)."), min_block), call. = FALSE)

  ## smallest block side caps the group order
  min_side <- min(vapply(blocks, function(b) min(length(b$rows), length(b$cols)), integer(1)))
  K_was_null <- is.null(K)
  K <- .default_K(K, min_side)

  ## The SMALLEST selected block is the binding constraint on resolution
  ## (audit F3.2): under the default K = min_side - 1, a small block caps the
  ## attainable p-value at 1/min_side even though it adds information. Say so
  ## explicitly when the cap makes rejection at this alpha impossible; the
  ## generic engine note states the arithmetic, this one names the cause and
  ## the remedy.
  res_note <- character(0)
  if (K_was_null && is.numeric(alpha) && length(alpha) == 1L &&
      is.finite(alpha) && alpha > 0 && alpha < 1 && 1 / (K + 1) > alpha) {
    res_note <- sprintf(paste0(
      "The smallest selected block (min side %d) caps the group order: ",
      "K = %d, so no rejection is attainable at alpha = %.3g (smallest ",
      "p-value 1/%d = %.3g). Testing at this level needs a fully observed ",
      "block with both sides >= %d; raise `min_block` to stop small blocks ",
      "from setting K."),
      min_side, K, alpha, K + 1L, 1 / (K + 1), ceiling(1 / alpha))
  }

  ## --- restrict data to cells inside the selected blocks ---------------------
  ## For each retained cell we record which block it belongs to and its position
  ## *within* that block (local row/col index), so the permutation builder can
  ## permute block-locally and translate back to global cluster ids.
  keep <- logical(N)                 # is this observation inside some block?
  blk_of <- integer(N)               # block index q for each observation
  lrow <- integer(N); lcol <- integer(N)   # 1-based position within block b$rows / b$cols
  for (q in seq_along(blocks)) {
    b <- blocks[[q]]
    in_rows <- ri %in% b$rows
    in_cols <- ci %in% b$cols
    sel <- in_rows & in_cols         # cells of this fully observed block
    keep <- keep | sel
    blk_of[sel] <- q
    lrow[sel] <- match(ri[sel], b$rows)
    lcol[sel] <- match(ci[sel], b$cols)
  }
  ## Subset every per-observation vector to the retained cells (suffix _k).
  idx <- which(keep)
  yk <- y[idx]; Dk <- D[idx, , drop = FALSE]; Xk <- X[idx, , drop = FALSE]
  ri_k <- ri[idx]; ci_k <- ci[idx]                      # global cluster ids, retained
  blk_k <- blk_of[idx]; lrow_k <- lrow[idx]; lcol_k <- lcol[idx]   # block + local indices
  code_k <- .cell_code(cbind(ri_k, ci_k))   # global cell code on retained cells
  Nk <- length(idx)                  # number of cells actually used

  sizes <- vapply(blocks, function(b) c(length(b$rows), length(b$cols)),
                  integer(2))
  note <- c(
    sprintf("Used %d of %d observed cells (%.1f%%) across %d fully observed block(s).",
            Nk, N, 100 * Nk / N, length(blocks)),
    sprintf("Block sizes (rows x cols): %s.",
            paste(sprintf("%dx%d", sizes[1, ], sizes[2, ]), collapse = ", ")),
    res_note
  )

  ## --- per-rep observation-permutation builder (block-diagonal) --------------
  perm_builder <- function(rep_seed)
    .build_obs_perms_blocks(rep_seed, K, blocks,
                            ri = ri_k, ci = ci_k, blk = blk_k,
                            lrow = lrow_k, lcol = lcol_k)

  res <- .ipt_engine(yk, Dk, Xk, perm_builder, K = K, n_reps = n_reps,
                     seed = seed, alpha = alpha, conf_int = conf_int,
                     beta_null = beta_null, grid = grid,
                     type = "missing (bicliques)", d_names = d_names,
                     n_clusters = c(row = n_row, col = n_col), call = cl,
                     n_cores = n_cores)
  res$note <- c(note, res$note)
  res$n_blocks <- length(blocks)
  res$cells_used <- Nk
  res$cells_total <- N
  res
}

#' Block-diagonal observation permutations for the missing-data design
#'
#' Permutes rows and columns independently within each fully observed block
#' (the blocks are disjoint in both dimensions, so the joint map is a genuine
#' relabelling), then expresses the result as observation gather-vectors over
#' the retained cells. Extracted from \code{mwperm_missing()}'s inline closure
#' so the group-closure tests exercise the SAME code the fit runs (audit
#' F2.6); behaviour-identical, same seed offsets.
#'
#' @param rep_seed rep-level seed (or NULL); block q draws its row/col groups
#'   at \code{.sub_seed(rep_seed, 4q - 1)} / \code{.sub_seed(rep_seed, 4q)}.
#' @param K group order minus one (common across blocks and dimensions).
#' @param blocks list of blocks as returned by \code{find_bicliques()}.
#' @param ri,ci global row/col cluster ids of the retained cells.
#' @param blk block index of each retained cell.
#' @param lrow,lcol 1-based position of each retained cell within its block's
#'   \code{rows}/\code{cols} vectors.
#' @return list of K+1 integer gather-vectors over the retained cells.
#' @keywords internal
#' @noRd
.build_obs_perms_blocks <- function(rep_seed, K, blocks, ri, ci, blk,
                                    lrow, lcol) {
  ## Per block, a row group and a col group of common order K+1. The 4*q
  ## offsets keep every block's two seeds distinct within a rep.
  rowG <- vector("list", length(blocks))
  colG <- vector("list", length(blocks))
  for (q in seq_along(blocks)) {
    rowG[[q]] <- build_perm_set(length(blocks[[q]]$rows), K,
                                seed = .sub_seed(rep_seed, 4L * q - 1L))
    colG[[q]] <- build_perm_set(length(blocks[[q]]$cols), K,
                                seed = .sub_seed(rep_seed, 4L * q))
  }
  code <- .cell_code(cbind(ri, ci))   # global cell code on retained cells
  Kp1 <- K + 1L
  ops <- vector("list", Kp1)          # one observation gather-vector per element
  for (k in seq_len(Kp1)) {
    ## Compute the target global (row, col) of every retained cell under the
    ## k-th element, block by block (unchanged where the block has no cells).
    tg_row <- ri; tg_col <- ci
    for (q in seq_along(blocks)) {
      sel <- blk == q                # retained cells in block q
      if (!any(sel)) next
      rimg <- rowG[[q]][[k]]; cimg <- colG[[q]][[k]]   # local image vectors
      tg_row[sel] <- blocks[[q]]$rows[rimg[lrow[sel]]]   # local -> permuted -> global
      tg_col[sel] <- blocks[[q]]$cols[cimg[lcol[sel]]]
    }
    ## Map each target cell back to its observation index among retained cells.
    tgt_code <- .cell_code(cbind(tg_row, tg_col))
    g <- match(tgt_code, code)
    if (anyNA(g))
      stop("Internal error: block permutation left the observed set.",
           call. = FALSE)   # nocov
    ops[[k]] <- g
  }
  ops
}


#' Disjoint fully observed bicliques (greedy or exact)
#'
#' Partitions (part of) a sparse two-way layout into blocks of row clusters
#' \eqn{I_q} and column clusters \eqn{J_q} such that every cell of
#' \eqn{I_q \times J_q} is observed and the blocks share no row or column
#' cluster. Blocks are peeled off one at a time: the largest fully observed
#' biclique among the still-available rows and columns is found, removed, and
#' the search repeats on the remainder. Used by \code{\link{mwperm_missing}}.
#'
#' Two ways to find each block are offered through \code{method}:
#' \describe{
#'   \item{\code{"greedy"}}{(default) a fast seed-and-intersect heuristic that,
#'     starting from each of a few highest-degree rows, adds rows while keeping
#'     the common set of fully observed columns. Always valid but may return a
#'     sub-maximal block.}
#'   \item{\code{"exact"}}{a branch-and-bound search for the maximum-area fully
#'     observed biclique, with area-based pruning. Maximum-edge biclique is
#'     NP-hard, so a per-block node budget caps the work; if it is hit the
#'     search falls back to the greedy block for that step (still valid) and a
#'     warning is issued. Note what is maximised: each block \emph{in turn},
#'     not the total covered area of the returned partition -- the greedy
#'     peeling after each exact block can leave less for later blocks, so
#'     \code{"exact"} does not always cover more cells than \code{"greedy"}
#'     overall, and it costs orders of magnitude more time on large masks.
#'     Worthwhile when a single largest block (finest p-value resolution)
#'     matters more than total coverage.}
#' }
#'
#' @param row,col Integer (or factor-coercible) cluster ids of the observed
#'   cells; the two vectors have equal length, one entry per observed cell.
#' @param min_block Integer; blocks whose smaller side is below this are not
#'   returned. Defaults to 2; values below 2 are silently raised to 2 (a
#'   one-sided block cannot be permuted).
#' @param method Either \code{"greedy"} (default) or \code{"exact"}; see
#'   Details.
#' @param node_budget Integer node cap for the \code{"exact"} branch-and-bound
#'   per block (default 200000). Ignored when \code{method = "greedy"}.
#'
#' @return A list of blocks, each a list with integer components \code{rows} and
#'   \code{cols} giving the (original-coding) cluster ids in that biclique.
#'
#' @seealso \code{\link{mwperm_missing}}.
#' @examples
#' ## a 4x4 grid missing its diagonal
#' g <- expand.grid(i = 1:4, j = 1:4)
#' g <- g[g$i != g$j, ]
#' find_bicliques(g$i, g$j, min_block = 2)
#' find_bicliques(g$i, g$j, min_block = 2, method = "exact")
#' @export
find_bicliques <- function(row, col, min_block = 2L,
                           method = c("greedy", "exact"),
                           node_budget = 2e5L) {
  method <- match.arg(method)
  min_block <- max(2L, as.integer(min_block))
  rf <- factor(row); cf <- factor(col)
  ri <- as.integer(rf); ci <- as.integer(cf)
  nr <- nlevels(rf); nc <- nlevels(cf)            # number of row / col clusters
  ## Observed-cell incidence as a logical matrix A[row, col] (cluster counts are
  ## usually modest, so a dense matrix is fine). A[i, j] == TRUE iff cell (i,j) seen.
  A <- matrix(FALSE, nr, nc)
  A[cbind(ri, ci)] <- TRUE

  ## Original cluster labels, restored to numeric where possible so the returned
  ## blocks use the user's own coding rather than factor levels.
  row_lab <- levels(rf); col_lab <- levels(cf)
  restore <- function(lab) {
    num <- suppressWarnings(as.numeric(lab))
    if (!anyNA(num)) num else lab
  }
  row_lab <- restore(row_lab); col_lab <- restore(col_lab)

  ## Peel off disjoint blocks one at a time: rows/cols already claimed by a block
  ## become unavailable for later blocks.
  avail_row <- rep(TRUE, nr); avail_col <- rep(TRUE, nc)
  blocks <- list()
  budget_hit <- FALSE                  # did the exact search ever exhaust its budget?

  repeat {
    rrows <- which(avail_row); rcols <- which(avail_col)   # still-available clusters
    if (length(rrows) < min_block || length(rcols) < min_block) break
    Asub <- A[rrows, rcols, drop = FALSE]   # incidence among available clusters
    if (method == "exact") {
      blk <- .max_biclique_exact(Asub, node_budget = node_budget)
      if (!isTRUE(blk$exact)) {        # budget hit: fall back to the greedy block
        budget_hit <- TRUE
        gre <- .grow_biclique(Asub)               # valid fallback for this block
        if (length(gre$rows) * length(gre$cols) > blk$area) blk <- gre
      }
    } else {
      blk <- .grow_biclique(Asub)         # local indices into rrows/rcols
    }
    if (length(blk$rows) < min_block || length(blk$cols) < min_block) break

    ## Map the block's local indices back to global clusters, store, and retire them.
    gr <- rrows[blk$rows]; gc <- rcols[blk$cols]
    blocks[[length(blocks) + 1L]] <- list(rows = row_lab[gr], cols = col_lab[gc])
    avail_row[gr] <- FALSE; avail_col[gc] <- FALSE
  }
  if (budget_hit)
    warning("Exact biclique search hit its node budget on at least one block; ",
            "used the greedy block there. Raise `node_budget` or use ",
            "method = \"greedy\".", call. = FALSE)
  blocks
}

#' Grow one large all-ones submatrix from a logical matrix (greedy)
#' @param A logical matrix.
#' @return list(rows, cols) of local indices forming an all-ones submatrix.
#' @keywords internal
#' @noRd
.grow_biclique <- function(A) {
  nr <- nrow(A); nc <- ncol(A)
  if (nr == 0L || nc == 0L) return(list(rows = integer(0), cols = integer(0)))
  deg <- rowSums(A)                        # per-row number of observed columns
  ord <- order(deg, decreasing = TRUE)     # process the densest rows first
  best <- list(rows = integer(0), cols = integer(0), area = 0)  # best block so far

  ## Try a few high-degree seed rows; keep the largest-area block found.
  n_seed <- min(nr, 8L)
  for (s in seq_len(n_seed)) {
    seed_row <- ord[s]
    cols <- which(A[seed_row, ])           # columns observed for the seed row
    if (length(cols) == 0L) next
    rows <- seed_row                       # current row set (grows below)
    area <- length(rows) * length(cols)    # current all-ones area
    ## Consider the remaining rows in degree order; adding a row shrinks the
    ## common column set to those it also observes. Keep a row only if the
    ## resulting block is at least as large (greedy area maximisation).
    for (r in ord) {
      if (r == seed_row) next
      new_cols <- cols[A[r, cols]]          # cols still fully observed if r added
      new_area <- (length(rows) + 1L) * length(new_cols)
      if (new_area >= area && length(new_cols) >= 1L) {
        rows <- c(rows, r); cols <- new_cols; area <- new_area
      }
    }
    if (area > best$area) best <- list(rows = rows, cols = cols, area = area)
  }
  list(rows = sort(best$rows), cols = sort(best$cols))
}

#' Maximum-area all-ones submatrix by branch-and-bound (exact)
#'
#' Finds the row set R maximising \code{|R| * |cols(R)|}, where
#' \code{cols(R)} is the set of columns observed for every row in R (so the
#' selected submatrix is all ones by construction). Rows are processed in
#' decreasing support order and the branch is pruned whenever the optimistic
#' bound \code{(rows so far + rows remaining) * current common columns} cannot
#' beat the incumbent. A node budget bounds the work; if exhausted the best
#' block found so far is returned with \code{exact = FALSE}.
#'
#' @param A logical matrix.
#' @param node_budget integer cap on the number of search nodes.
#' @return list(rows, cols, area, exact); rows/cols are local indices.
#' @keywords internal
#' @noRd
.max_biclique_exact <- function(A, node_budget = 2e5L) {
  nr <- nrow(A); nc <- ncol(A)
  if (nr == 0L || nc == 0L)
    return(list(rows = integer(0), cols = integer(0), area = 0, exact = TRUE))
  supp <- lapply(seq_len(nr), function(r) which(A[r, ]))   # column support per row
  ord  <- order(lengths(supp), decreasing = TRUE)          # process dense rows first
  best <- list(rows = integer(0), cols = integer(0), area = 0)  # incumbent best block
  nodes <- 0L                          # search nodes visited (vs. node_budget)
  exact <- TRUE                        # cleared to FALSE if the budget is hit

  ## Depth-first branch-and-bound over rows in `ord`.
  ##   pos     : index into `ord` of the row being decided
  ##   chosen  : row indices included so far
  ##   curcols : columns observed by *every* chosen row (the block's columns)
  rec <- function(pos, chosen, curcols) {
    if (nodes >= node_budget) { exact <<- FALSE; return(invisible(NULL)) }
    nodes <<- nodes + 1L
    rem <- length(ord) - pos + 1L      # rows not yet decided
    if (length(chosen) > 0L) {
      ncols <- length(curcols)
      if (ncols == 0L) return(invisible(NULL))           # no common column: dead end
      area <- length(chosen) * ncols
      if (area > best$area) best <<- list(rows = chosen, cols = curcols, area = area)
      ## Optimistic bound: even adding every remaining row keeps <= ncols columns,
      ## so if that cannot beat the incumbent, prune this branch.
      if ((length(chosen) + rem) * ncols <= best$area) return(invisible(NULL))
    }
    if (pos > length(ord)) return(invisible(NULL))
    r <- ord[pos]
    ## Branch 1: include row r (intersect its support with the running columns).
    newcols <- if (length(chosen) == 0L) supp[[r]] else curcols[curcols %in% supp[[r]]]
    rec(pos + 1L, c(chosen, r), newcols)
    if (!exact) return(invisible(NULL))   # budget hit deeper down: stop early
    ## Branch 2: skip row r.
    rec(pos + 1L, chosen, curcols)
  }
  rec(1L, integer(0), integer(0))
  list(rows = sort(best$rows), cols = sort(best$cols), area = best$area,
       exact = exact)
}
