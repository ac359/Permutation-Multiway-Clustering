## Sign-flip (Rademacher) invariant test for dyadic regression: IPT-Het.
##
## The permutation designs rest on exchangeability of the error array under
## relabelling of the clusters (Assumption 1 of Guo, Toulis and Wang 2026).
## That fails under heteroskedasticity, because relabelling clusters
## relabels the variance pattern. Section 2 of the paper notes that the
## partialling-out + minorization argument is not specific to permutations
## and goes through under ANY invariance group G. This file cashes that in
## with a sign-flip group, which needs a different, non-nested assumption:
##
##     (eps_ij) =d (s_i t_j eps_ij) | X, D     for row signs s, column signs t.
##
## A sign flip changes no variance, so arbitrary heteroskedasticity -- in i,
## in j, in the covariates -- is fine; the price is symmetry of the errors
## about zero, which exchangeability never asked for.
##
## The group. Fix n_flip. Each row cluster and each column cluster is
## assigned independently and uniformly to one of n_flip flip groups (g1, g2).
## A sign vector s in {-1, +1}^n_flip acts on observation (i, j) by the sign
## s[g1(i)] * s[g2(j)]; write S_s for that diagonal +/-1 matrix. The set
## {S_s} is a group (S_s S_s' = S_{s * s'}, with * the coordinate-wise
## product), so Theorem 1 applies with G = {S_s}, and the engine's Procedure 1
## core runs on it unchanged: element k residualizes on [X | S_k X] and the
## statistics a_k(b) = ||D' V_k V_k' (y - D b)||, b_k(b) = ||D' V_k V_k' S_k
## (y - D b)|| are affine in b exactly as in the permutation case, so the
## whole confidence-set path is inherited.
##
## The kernel, and why deduplication is exact. Because the sign enters as a
## PRODUCT of a row sign and a column sign, s and -s induce the identical
## S_s: the action has a kernel {s, -s} of order 2, so the 2^n_flip sign
## vectors produce only 2^(n_flip - 1) distinct transformations, each of them
## twice. Enumerating one representative per coset -- coordinate 1 fixed to
## +1 -- gives every distinct element exactly once. That leaves the p-value
## bit-identical: with the full enumeration every a_k, b_k appears twice, so
## min_j a_j is unchanged and both the count and the group order double.
## The effective group order is therefore 2^(n_flip - 1), the smallest
## attainable p-value 1/2^(n_flip - 1), and a 95% confidence set needs
## n_flip >= 6 (order 32 >= 20).
##
## The kernel is EXACTLY of order 2 only if every flip group is used by at
## least one row or column cluster (build_flip_set() resamples until it is):
## an unused group's coordinate is free, the kernel grows to order 4, and the
## 2^(n_flip - 1) representatives would then contain duplicates -- the test
## would still be valid, but the reported resolution 1/2^(n_flip - 1) would
## be a lie, because the true number of distinct elements is smaller.
##
## The identity is EXCLUDED from min_j a_j, as Equation (10) of the paper
## writes it (1 <= j <= K over the non-identity elements) and as the author
## confirmed for this design. The engine already does this for every design.
##
## Cost. One residual projection per non-identity element per repetition,
## i.e. 2^(n_flip - 1) - 1 of them: 31 at n_flip = 6, 127 at the default 8,
## 524,287 at the cap of 20. Exponential in n_flip where the permutation
## designs are linear in K; hence the fixed default and the cap
## (.default_n_flip() in engine.R).

#' Construct a random sign-flip group (the sign-flip analogue of Algorithm 1)
#'
#' Builds the `2^(n_flip - 1)` distinct sign-flip transformations of an
#' `n_row` by `n_col` array that [mwperm_dyadic_het()] uses as its invariance
#' group. Each row cluster and each column cluster is assigned, independently
#' and uniformly at random, to one of `n_flip` *flip groups*; a sign vector
#' `s` in `{-1, +1}^n_flip` then multiplies observation `(i, j)` by
#' `s[g1(i)] * s[g2(j)]`, the product of its row group's sign and its column
#' group's sign. The returned set is closed under composition (composition is
#' the coordinate-wise product of sign vectors), contains the identity
#' exactly once, and has exactly `2^(n_flip - 1)` elements.
#'
#' **Why `2^(n_flip - 1)` and not `2^n_flip`.** Because the sign enters as a
#' *product* of a row sign and a column sign, `s` and `-s` induce the identical
#' transformation: the action has a kernel of order 2, so the `2^n_flip` sign
#' vectors produce only `2^(n_flip - 1)` distinct transformations, each of
#' them twice. This function enumerates one representative per pair
#' (coordinate 1 fixed to `+1`), which is exact: the p-value is bit-identical
#' to what the full enumeration would give, at half the cost. The effective
#' group order is therefore `2^(n_flip - 1)` and the smallest attainable
#' p-value `1 / 2^(n_flip - 1)`, so a 95% confidence set needs `n_flip >= 6`.
#'
#' **Every flip group must be used.** If some group were assigned to no row
#' cluster and no column cluster its sign would be free, the kernel would be
#' larger than 2, and the `2^(n_flip - 1)` representatives would contain
#' duplicates -- a valid test whose *reported* resolution would be wrong. The
#' assignment is therefore resampled until the union of the row and column
#' images covers all `n_flip` groups. That is rare at small `n_flip` (about
#' 5% of draws at `n_flip = 10` with 25 clusters per side) but it must be
#' guaranteed, not left to chance.
#'
#' The randomness in the assignment is what makes the resulting test a
#' *random* invariant test, exactly as the relabelling does in
#' [build_perm_set()]: different seeds give slightly different p-values, and
#' the `n_reps` argument of [mwperm_dyadic_het()] aggregates several draws.
#' The same reproducibility caveats apply (same `RNGkind()`; see
#' [build_perm_set()]).
#'
#' @param n_row,n_col Integer cluster counts along the two dimensions.
#' @param n_flip Integer, the number of flip groups, at least 2 and at most
#'   `n_row + n_col` (so that every group can be reached). Group order is
#'   `2^(n_flip - 1)`, and the cost of the test is one residual projection
#'   per element, so keep it moderate: [mwperm_dyadic_het()] defaults to 8 and
#'   refuses values above 20.
#' @param seed Optional integer seed for the random assignment. If `NULL`,
#'   the current RNG state is used. A seeded call leaves the caller's RNG
#'   stream exactly as it found it.
#'
#' @return A list of length `2^(n_flip - 1)`. Element `k` is a list with two
#'   numeric `+1/-1` vectors, `row` (length `n_row`) and `col` (length
#'   `n_col`): the transformation multiplies observation `(i, j)` by `row[i] *
#'   col[j]`. The first element is the identity (all `+1`). Attributes:
#'   `"n_flip"`, `"group_order"` (`2^(n_flip - 1)`), `"row_groups"` and
#'   `"col_groups"` (the two integer assignments, length `n_row` and
#'   `n_col`), and `"signs"` (the `2^(n_flip - 1)` by `n_flip` matrix of coset
#'   representatives, one row per element, row 1 all `+1`).
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data. arXiv:2601.08610. Section 2
#'   (the argument holds for any invariance group) and Procedure 1.
#'
#' @seealso [mwperm_dyadic_het()], the test built on this group;
#'   [build_perm_set()] for the permutation group the other designs use.
#' @examples
#' F <- build_flip_set(n_row = 6, n_col = 5, n_flip = 4, seed = 1)
#' length(F)                  # 2^(4 - 1) = 8 distinct sign patterns
#' F[[1]]                     # identity: all +1
#' attr(F, "row_groups")      # which flip group each row cluster belongs to
#' attr(F, "group_order")
#' ## the sign applied to cell (i, j) by element 3
#' outer(F[[3]]$row, F[[3]]$col)
#' @export
build_flip_set <- function(n_row, n_col, n_flip, seed = NULL) {
  n_row  <- as.integer(n_row)
  n_col  <- as.integer(n_col)
  n_flip <- as.integer(n_flip)
  if (length(n_row) != 1L || is.na(n_row) || n_row < 1L ||
      length(n_col) != 1L || is.na(n_col) || n_col < 1L)
    stop("`n_row` and `n_col` must be single integers >= 1.", call. = FALSE)
  if (length(n_flip) != 1L || is.na(n_flip) || n_flip < 2L)
    stop("`n_flip` must be a single integer >= 2.", call. = FALSE)
  if (n_flip > n_row + n_col)
    stop(sprintf(paste0("`n_flip` = %d exceeds the number of clusters (%d ",
                        "rows + %d columns): every flip group must be used ",
                        "by at least one cluster. Use n_flip <= %d."),
                 n_flip, n_row, n_col, n_row + n_col), call. = FALSE)
  if (n_flip > 30L)                    # 2^30 elements; the front end caps at 20
    stop("`n_flip` must be at most 30.", call. = FALSE)

  ## Reproducible assignment without disturbing the caller's RNG stream (the
  ## same discipline as build_perm_set(): snapshot, reseed, restore on exit).
  if (!is.null(seed)) {
    old <- .save_seed()
    on.exit(.restore_seed(old), add = TRUE)
    set.seed(seed)
  }

  ## Row and column assignments, redrawn until every flip group is used by at
  ## least one cluster (see the header: an unused group enlarges the kernel
  ## and the 2^(n_flip - 1) representatives would no longer be distinct).
  ## The draw is `sample.int(n_flip, n, replace = TRUE)`, i.e. the same draw
  ## as `sample(n_flip, n, replace = TRUE)` in the reference script.
  for (attempt in seq_len(10000L)) {
    g1 <- sample.int(n_flip, n_row, replace = TRUE)
    g2 <- sample.int(n_flip, n_col, replace = TRUE)
    if (length(unique(c(g1, g2))) == n_flip) break
    if (attempt == 10000L)
      stop(sprintf(paste0("Could not find an assignment using all %d flip ",
                          "groups in 10000 draws (%d rows, %d columns); ",
                          "use a smaller `n_flip`."),
                   n_flip, n_row, n_col), call. = FALSE)
  }

  ## One representative per coset {s, -s}: coordinate 1 is pinned to +1 and
  ## the remaining n_flip - 1 coordinates run over {+1, -1}. expand.grid()
  ## varies its first factor fastest and starts at the first level, so with
  ## levels c(1, -1) row 1 is all +1: the identity comes first.
  S <- cbind(1, as.matrix(expand.grid(rep(list(c(1, -1)), n_flip - 1L),
                                      KEEP.OUT.ATTRS = FALSE)))
  dimnames(S) <- NULL
  storage.mode(S) <- "double"
  order <- nrow(S)                     # 2^(n_flip - 1)

  ## Element k acts on row cluster i by S[k, g1[i]] and on column cluster j by
  ## S[k, g2[j]]; the observation-level sign is their product.
  flips <- vector("list", order)
  for (k in seq_len(order))
    flips[[k]] <- list(row = S[k, g1], col = S[k, g2])
  attr(flips, "n_flip")      <- n_flip
  attr(flips, "group_order") <- order
  attr(flips, "row_groups")  <- g1
  attr(flips, "col_groups")  <- g2
  attr(flips, "signs")       <- S
  flips
}

#' Build observation-level sign-flip operators from a flip set
#'
#' The sign-flip counterpart of `.build_obs_perms()`: composes the per-cluster
#' row and column signs of a [build_flip_set()] element into one `+1/-1`
#' vector over the observations and wraps it as the signed gather `list(g =
#' NULL, s = <signs>)` that `.apply_op()` and hence `.ipt_prepare()` accept.
#' There is no gather (`g = NULL`): a sign flip moves no observation, which is
#' also why it needs no completeness of the array and no bijection check --
#' the analogue of `.assert_bijection()` here is that every `s` has length N
#' and entries in `{-1, +1}`, asserted below.
#'
#' @param coords integer matrix N x 2 of dense (row, col) cluster ids.
#' @param flips a flip set from [build_flip_set()].
#' @param design short label for the calling design, used in error messages.
#' @return list of length `2^(n_flip - 1)` of signed gathers, identity first.
#' @keywords internal
#' @noRd
.build_obs_flips <- function(coords, flips, design = "this design") {
  coords <- as.matrix(coords)
  N <- nrow(coords)
  ri <- coords[, 1L]
  ci <- coords[, 2L]
  ops <- vector("list", length(flips))
  for (k in seq_along(flips)) {
    s <- flips[[k]]$row[ri] * flips[[k]]$col[ci]   # sign of observation (i, j)
    if (length(s) != N || anyNA(s) || !all(s == 1 | s == -1))
      stop(sprintf(paste0("Internal error: element %d of the %s sign-flip ",
                          "group is not a +1/-1 vector over the %d ",
                          "observations. Please report this with a ",
                          "reproducible example."), k, design, N),
           call. = FALSE)
    ops[[k]] <- list(g = NULL, s = s)
  }
  ops
}

#' Sign-flip invariant test for dyadic regression (heteroskedasticity-robust)
#'
#' Finite-sample valid test of H0: beta = b in the dyadic regression model
#' \preformatted{  y_ij = x_ij' gamma + d_ij' beta + eps_ij}
#'
#' with a single observation per cell, under an invariance assumption that
#' tolerates **arbitrary heteroskedasticity**: the errors must be
#' symmetrically distributed under joint row-and-column sign changes,
#' \preformatted{  (eps_ij) =d (s_i t_j eps_ij) | X, D   for all row signs s and column signs t.}
#'
#' This is the invariant test of Guo, Toulis and Wang (2026) run with a
#' *sign-flip* group in place of the permutation group -- Section 2 of the
#' paper notes that the partialling-out and minorization argument holds for
#' any invariance group -- and is the test the authors call IPT-Het. Same
#' model, same null, same statistic, same confidence set by test inversion;
#' only the group changes, and with it the assumption.
#'
#' **When to use it instead of [mwperm_dyadic()].** The two assumptions are
#' complementary and neither contains the other:
#' - [mwperm_dyadic()] needs the error array to be *exchangeable* under
#'   relabelling of the clusters, conditional on the covariates. That fails
#'   whenever the error variance depends on the cluster identity or on the
#'   covariates -- a gravity equation whose residual variance grows with
#'   distance or GDP, say -- because relabelling the clusters relabels the
#'   variance pattern. Symmetry is *not* required.
#' - `mwperm_dyadic_het()` needs the errors to be *symmetric about zero*
#'   under row and column sign changes. A sign flip changes no variance, so
#'   any pattern of heteroskedasticity is fine; skewed errors are not.
#'
#' Under covariate-dependent heteroskedasticity (25 clusters per side, a
#' gravity design with the residual variance increasing in the mean and in
#' the distance covariate) the permutation test's rejection rate under a true
#' null rose with the strength of the heteroskedasticity -- to about
#' 0.06-0.09 at a nominal 0.05 in this package's own 300-replication check,
#' and to about 0.13 in the authors' -- while this test stayed at or below
#' nominal (0.02 and 0.05 respectively). That gap is the reason the function
#' exists. It is **not a strict upgrade**: under homoskedastic, exchangeable
#' errors the sign-flip test is less powerful (about 0.86 against 0.96 at
#' beta = 0.15 in the same design here; 0.90 against 0.98 in the authors'
#' run), because its group is smaller and coarser than a full relabelling
#' group. Use it when you have reason to doubt exchangeability and are
#' prepared to assume symmetry; use [mwperm_dyadic()] otherwise.
#' Heteroskedasticity leaves no trace in the clustering structure, so
#' [mwperm()] never selects this test automatically: request it with
#' `design = "dyadic_het"`.
#'
#' **The group and its resolution.** Each row cluster and each column cluster
#' is assigned at random to one of `n_flip` flip groups, and a sign vector in
#' `{-1, +1}^n_flip` multiplies observation `(i, j)` by the product of its
#' row group's sign and its column group's sign (see [build_flip_set()]). A
#' sign vector and its negative induce the same transformation, so the group
#' has `2^(n_flip - 1)` distinct elements, not `2^n_flip`; the package
#' enumerates each exactly once. Per repetition the p-value therefore lives
#' on the grid `{1, ..., 2^(n_flip - 1)} / 2^(n_flip - 1)`, its smallest
#' value is `1 / 2^(n_flip - 1)`, and a `(1 - alpha)` confidence set needs
#' `2^(n_flip - 1) >= 1 / alpha`: at `alpha = 0.05` that is `n_flip >= 6`
#' (order 32). Below that the p-value is still exact but no set is
#' attainable, and the fit says so in a note. The cost is one residual
#' projection per non-identity element -- `2^(n_flip - 1) - 1` per
#' repetition, exponential in `n_flip` -- which is why the default is a fixed
#' `n_flip = 8` (order 128, floor 0.0078, 127 projections) rather than the
#' largest value the design allows, as `K` is for the permutation tests. At
#' that default with `n_reps = 10` the exact confidence set exceeds the
#' engine's candidate budget and the interval is found by bracketing and
#' bisection instead (the fit says so in a note; the end points agree to the
#' bisection tolerance); `n_reps <= 6`, or `n_flip <= 7`, keeps it exact.
#'
#' The aggregation over `n_reps` repetitions, the exact confidence set and
#' its closure convention, the `aggregate = "median2"` guarantee and the
#' numerical notes on the projection are all exactly as documented in
#' [mwperm_dyadic()]; the sign-flip element `S_k` simply takes the place of
#' the permuted copy in the stacked design `[X | S_k X]`.
#'
#' @inheritParams mwperm_dyadic
#' @param n_flip Number of flip groups. The group has order `2^(n_flip - 1)`,
#'   so the smallest attainable p-value is `1 / 2^(n_flip - 1)` and a 95%
#'   confidence set needs `n_flip >= 6`. Defaults to 8, capped at `min(n_row,
#'   n_col)` (every group must be reachable by the row assignment alone);
#'   values above 20 are refused, with the number of projections they would
#'   imply. Must satisfy `2 <= n_flip <= min(n_row, n_col)`.
#' @param n_reps Number of independent runs whose p-values are aggregated by
#'   `aggregate` (the median by default), as recommended for randomised
#'   tests; the confidence set inverts the same aggregated p-value. Defaults
#'   to 10: a single run's p-value depends on the random assignment of
#'   clusters to flip groups (a seed lottery), and the median of 10 runs
#'   stabilises it at roughly ten times the cost. See *Aggregation over
#'   repetitions* in [mwperm_dyadic()].
#' @param aggregate How the `n_reps` per-repetition p-values are combined
#'   into the reported p-value, and into the confidence set that inverts it.
#'   `"median"` (the default) is the median, as recommended in Remark 1 of
#'   Guo, Toulis and Wang (2026); `"median2"` is `min(1, 2 * median)`, a
#'   valid p-value at level `alpha` under arbitrary dependence across
#'   repetitions, at the cost of resolution: its floor is
#'   `2 / 2^(n_flip - 1)`, so a 95% set then needs `n_flip >= 7`. See
#'   [mwperm_dyadic()] for the full account.
#' @param n_cores Number of CPU cores (default 1 = serial). Parallelism is
#'   over the `n_reps` repetitions when several are run with a `seed`,
#'   otherwise over the `2^(n_flip - 1) - 1` per-element projections; either
#'   way the result is *identical* to the serial one. See [mwperm_dyadic()].
#'
#' @return An object of class `"mwperm"` exactly as [mwperm_dyadic()] returns
#'   it (see there for the provenance of every field), with `type` equal to
#'   `"dyadic (sign-flip / heteroskedasticity-robust)"`, `K` and `n_perm`
#'   equal to `2^(n_flip - 1) - 1` and `2^(n_flip - 1)` (the non-identity
#'   count and the group order, the roles those fields play for every
#'   design), and an extra field `n_flip`. All the methods -- [print.mwperm()],
#'   [summary.mwperm()], [confint.mwperm()], [plot.mwperm()] -- apply.
#'
#' @references Guo, W., Toulis, P. and Wang, Y. (2026). Permutation inference
#'   under multi-way clustering and missing data. arXiv:2601.08610.
#'
#' @seealso [mwperm_dyadic()] for the permutation test under exchangeability,
#'   the test to compare this one against; [build_flip_set()] for the group
#'   construction; [mwperm()] with `design = "dyadic_het"`.
#'
#' @examples
#' data(trade_dyadic)
#' ## the same call as the mwperm_dyadic() example, under sign symmetry
#' ## instead of exchangeability
#' fit <- with(trade_dyadic,
#'             mwperm_dyadic_het(y = log_trade, d = log_dist,
#'                               x = cbind(log_gdp_i, log_gdp_j),
#'                               row = importer, col = exporter,
#'                               n_flip = 6, n_reps = 3, seed = 1))
#' fit
#' confint(fit)
#' @export
mwperm_dyadic_het <- function(y, d, x = NULL, row, col, n_flip = NULL,
                              alpha = 0.05, beta_null = 0, conf_int = TRUE,
                              n_reps = 10L, seed = NULL, grid = NULL,
                              aggregate = c("median", "median2"),
                              n_cores = 1L) {
  cl <- match.call()                   # stored on the result for printing
  aggregate <- match.arg(aggregate)
  y <- .check_y(y)
  N <- length(y)                       # number of observations
  D <- as.matrix(d)
  d_names <- .coef_names(D, deparse(substitute(d)))  # coefficient label(s)
  X <- .make_X(x, N)                   # nuisance design with intercept
  .check_lengths(N, list(row = row, col = col))

  ri <- .dense_id(row, "row")
  ci <- .dense_id(col, "col")          # dense 1-based row/col cluster ids
  n_row <- max(ri)
  n_col <- max(ci)                     # number of row / col clusters
  if (anyDuplicated(cbind(ri, ci)))
    stop("Dyadic regression expects one observation per (row, col) cell. ",
         "For repeated observations use mwperm_layout() or mwperm_panel().",
         call. = FALSE)
  coords <- cbind(ri, ci)              # per-observation (row, col) coordinates
  ## Same contract as mwperm_dyadic(): one observation per cell of a complete
  ## array. (A sign flip moves no observation, so the construction itself
  ## would run on an incomplete array; the contract is kept identical to the
  ## permutation front end so the two are interchangeable in a call.)
  .require_complete_array(coords, c(row = n_row, col = n_col), N,
                          what = "Dyadic regression",
                          remedy = paste0("mwperm_missing(), which restricts ",
                                          "to fully observed blocks and ",
                                          "permutes within them"))

  n_flip <- .default_n_flip(n_flip, n_row, n_col)
  ## The engine's K is the non-identity count, whatever the construction:
  ## here 2^(n_flip - 1) - 1, so that its resolution guard sees the group
  ## order n_perm = 2^(n_flip - 1) and gates the confidence set on
  ## 1 / 2^(n_flip - 1).
  K_eff <- as.integer(2^(n_flip - 1L)) - 1L

  ## Per-rep group: one seeded assignment of clusters to flip groups (both
  ## dimensions from the one sub-seed, drawn in sequence), composed into
  ## observation-level sign vectors. Sub-seed offset 1, as the permutation
  ## front ends use for their first dimension; the two families of front end
  ## never share a seed scheme, so nothing existing is affected.
  perm_builder <- function(rep_seed) {
    F <- build_flip_set(n_row, n_col, n_flip, seed = .sub_seed(rep_seed, 1L))
    .build_obs_flips(coords, F, design = "dyadic (sign-flip)")
  }

  res <- .ipt_engine(y, D, X, perm_builder, K = K_eff, n_reps = n_reps,
                     seed = seed, alpha = alpha, conf_int = conf_int,
                     beta_null = beta_null, grid = grid,
                     type = "dyadic (sign-flip / heteroskedasticity-robust)",
                     d_names = d_names,
                     n_clusters = c(row = n_row, col = n_col), call = cl,
                     n_cores = n_cores, ci_agg = aggregate, group = "flip")
  res$n_flip <- n_flip                 # read by print() / confint()
  res
}
