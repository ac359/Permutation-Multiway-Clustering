# mwperm architecture

This file is for someone who knows the method paper and R, and has never
opened this code. It maps each step of the invariant permutation test to the
code, names the variable behind each paper symbol, and lists the conventions a
correct change must keep. It is listed in `.Rbuildignore`, so it is not in
the package tarball.

**Sources.**

* **GTW:** Guo, Toulis & Wang (2026), *Permutation inference under
  multi-way clustering and missing data*, arXiv:2601.08610. The sign-flip
  test follows the revised paper's Assumption 2 (double sign symmetry).
* **WWW:** Wen, Wang & Wang (2025), *Residual permutation test for
  regression coefficient testing*, Ann. Statist. 53(2).

`TESTING_PLAN.md` holds the traceability table (paper object -> function ->
tests) and the log of every place where code and paper differ.

## 1. The pipeline

```
mwperm() / mwperm_formula()                         R/unified.R, R/formula.R
  |  mwperm_check(): read the index structure, pick the design
  v
design worker: mwperm_dyadic | _panel | _threeway | _layout | _missing |
               _panel_missing | _irregular | _dyadic_het
  |  validate, build X (intercept + x [+ period dummies]), fix K,
  |  and define perm_builder(rep_seed) -> K + 1 group elements
  v
.ipt_engine()                                                   R/engine.R
  |  for r in 1..n_reps (seed + r - 1):
  |    permutation construction   perm_builder(rep_seed)
  |      build_perm_set()      Algorithm 1, one dimension     R/perm_set.R
  |      .build_obs_perms()    Eq. (9), lifted to rows        R/core.R
  |      (.build_obs_perms_layout, .build_obs_perms_blocks,
  |       build_flip_set + .build_obs_flips for the variants)
  |    projection engine
  |      .ipt_prepare()        Procedure 1, step 1            R/core.R
  |      .ipt_eval()           Procedure 1, step 2 = Eq. (10) R/core.R
  |  median aggregation
  |    .agg_pvals()            Remark 1                       R/engine.R
  |  test inversion
  |    .invert_ci()            step 3, d = 1: .ci_breakpoints + .exact_ci_set
  |    .invert_region()        step 3, d > 1 (grid)
  v
"mwperm" object -> print / summary / confint / coef / nobs / plot
                                           R/methods.R, R/formula.R, R/plot.R
```

A fit touches the data twice. Once per group element, in `.ipt_prepare()`: the
residual projection, which is the expensive part. Once per candidate null, in
`.ipt_eval()` / `.pval_matrix()`: arithmetic on cached d x d quantities, since
the statistic is affine in the null `b`. So the confidence set costs no extra
matrix factorization.

## 2. Paper notation -> code

| Paper (GTW / WWW) | Code | Where |
|---|---|---|
| y, D, X in model (7); X includes the intercept | `y`, `D` (N x d matrix), `X` (N x p) | every front end; `.make_X()` |
| stacking cell (i, j) in row (i-1) n + j | not used: cells are keyed by `.cell_code()` (mixed radix, first coordinate least significant); the statistic does not depend on row order | `R/core.R` |
| n row clusters, m / n / ell dimension sizes | `n_row`, `n_col`, `m`, `n`, `ell`; dense ids `ri`, `ci` | front ends; `.dense_id()` |
| K; group order K + 1 | `K`; `Kp1`, `n_perm` | `.default_K()`, engine |
| Algorithm 1: pi, psi~_k, psi_k | `pi_vec`, `psit`, `perms[[k + 1]]` | `build_perm_set()` |
| (pi_k, sigma_k) of Eq. (9) | `Grow[[k + 1]]`, `Gcol[[k + 1]]` | `mwperm_dyadic()` perm_builder |
| y_{pi_k,sigma_k}, X_{pi_k,sigma_k} | `y[g]`, `X[g, ]` with gather vector `g = obs_perms[[k + 1]]`, applied by `.apply_op()` | `.build_obs_perms()`, `.apply_op()` |
| S_k (sign-flip element) | `list(g = NULL, s = <+/-1>)` | `.build_obs_flips()` |
| V_k V_k' (projector onto col([X, X_k])^perp); WWW: V~_k V~_k' | never formed; `Dr = D - [X, X_k] cf` is V_k V_k' D | `.ipt_prepare()` |
| a_k = norm(D' V_k V_k' y) | `u[, k]` at b = 0; `a` in `.ipt_eval()` | `R/core.R` |
| b_k = norm(D' V_k V_k' y_k) | `v[, k]` at b = 0; `b` in `.ipt_eval()` | `R/core.R` |
| slopes in the null b | `M[, , k]` = D'V V'D, `W[, , k]` = D'V V'D_k | `R/core.R` |
| min_{1 <= j <= K} a_j | `amin` | `.ipt_eval()`, `.pval_matrix()` |
| pval, Eq. (10) | `(1 + sum(b >= amin)) / Kp1` | `.ipt_eval()` |
| pval_r for repetition r; the median (Remark 1) | `pvalues_rep`; `pvalue` | `.ipt_engine()`, `.agg_pvals()` |
| CI = {b : pval(b) > alpha} (step 3) | `conf_set` (components), `conf_int` (hull), `conf_region` | `.invert_ci()`, `.invert_region()` |
| mask M, blocks F_M = {I_q x J_q} (Definition 1) | `blocks[[q]]$rows`, `$cols` | `find_bicliques()` |
| block-local pi_{k,q}, sigma_{k,q} | `rowG[[q]][[k]]`, `colG[[q]][[k]]` | `.build_obs_perms_blocks()` |
| replicate index l, cell size ell_ij, threshold L0 | `widx`/`slot`, `ell`, `L0` | `R/layout.R`, `R/irregular.R` |
| period t (held fixed under InvB) | the NULL group in `.build_obs_perms()`; `slot` in the block builder | `R/panel.R`, `R/panel_missing.R` |

## 3. File -> paper map

| File | Implements | Pipeline stage |
|---|---|---|
| `R/perm_set.R` | GTW Algorithm 1, Proposition 2; RNG helpers | permutation construction |
| `R/core.R` | Eq. (9) at the observation level; Procedure 1 steps 1-2, Eq. (10) | permutation construction, projection engine |
| `R/engine.R` | the repetition loop; Remark 1; Procedure 1 step 3 (exact set, grid, bisection, joint region); Theorem 1's p < N/2; seeds | projection engine, aggregation, inversion |
| `R/dyadic.R` | Section 3: Procedure 1 under Assumption 1 (InvA) | design worker |
| `R/threeway.R` | Section 6.1 (InvA, three indices) | design worker |
| `R/panel.R` | Section 6.2 (InvB) | design worker |
| `R/layout.R` | Section 6.3; the L0 threshold of Section 6.4 | design worker, construction |
| `R/missing.R` | Section 5: Assumption 4, Definition 1, Procedure 2, Theorem 4; Appendix A: Algorithm 2 | design worker, construction |
| `R/panel_missing.R` | Section 6.2 combined with Procedure 2 (the case GTW Section 9 leaves open) | design worker |
| `R/irregular.R` | Section 6.4, steps (i)-(ii) | design worker |
| `R/signflip.R` | Procedure 1 under the revised Assumption 2 (double sign symmetry) | design worker, construction |
| `R/unified.R` | the dispatch rule (which procedure the data support) | dispatch |
| `R/formula.R` | argument assembly for model (7); `coef()`, `nobs()` | entry, S3 methods |
| `R/methods.R`, `R/plot.R` | reporting only | S3 methods |
| `R/data.R`, `R/mwperm-package.R` | documentation; `NAMESPACE` imports | -- |

Every function in `R/` has a roxygen block whose `@details` line cites the
step it implements, or says "not a paper step". Every file opens with a header
giving its purpose, its paper sections and its pipeline stage.

## 4. Conventions a correct change must keep

1. **Validity comes from the group.** Each `perm_builder` returns a list of
   K + 1 elements, the identity first, and the list must be closed under
   composition (Theorem 1). `build_perm_set()` guarantees closure: its
   elements are powers of one generator. Every observation-level builder
   ends in `.assert_bijection()`.
2. **The group-element contract is `.apply_op()`'s.** An element is either a
   bare gather vector, applied as `M[g, ]`, or a signed gather
   `list(g, s)`. `.ipt_prepare()` may rely only on properties of orthogonal
   row actions.
3. **Rank, not 2p.** `[X | X_k]` is rank-deficient by construction: the
   intercept and panel period dummies map to themselves. The projection
   uses the numerical rank (eigenvalue cut 1e-14 of the Gram matrix).
   Do not "fix" this into N - 2p.
4. **Eq. (10) exactly.** The code keeps all of these:
   * the minimum runs over the K non-identity elements;
   * the indicator is `<=`, so ties count;
   * p-values live on `{1, ..., K+1}/(K+1)`;
   * the decision is `p <= alpha`.
5. **One aggregation rule** (`.agg_pvals()`) serves both the reported
   p-value and every confidence-set path, so the test and the set cannot
   disagree. The exact set's components are reported closed. This is the
   conservative side, and it is documented.
6. **Frozen seed scheme.** Repetition r uses `seed + r - 1`. Dimension, cell
   or block j uses `rep_seed * stride + j` (`.sub_seed()`); block q uses
   offsets 4q - 1 and 4q. Seeded draws save and restore `.Random.seed`.
   Changing any of this moves every seeded number in the golden baseline
   and the paper.
7. **Parallel equals serial.** Only seeded or RNG-free loops run in
   parallel. With `seed = NULL`, the repetition loop stays serial.

## 5. Adding a design

Write a new front end that validates its input and builds `X`. Then define
`perm_builder(rep_seed)`, returning K + 1 gather vectors (or signed gathers)
that form a group under the new invariance, and call `.ipt_engine()`. Do not
modify the engine. `R/panel_missing.R` is the worked example: it reuses
`find_bicliques()` and `.build_obs_perms_blocks(slot = )` and adds no engine
code.

A new design also ships with:

* a reference check in `tests/testthat/test-designs.R` (its p-values against
  the naive Procedure 1 of `helper-reference.R`, with the group rebuilt from
  the seed);
* a Monte Carlo size check in `tests/testthat/test-montecarlo.R`;
* a row in the traceability table of `TESTING_PLAN.md`.

## 6. Tests

`tests/README.md` describes the two suites. `tests/test-*.R` and
`tests/lower-level-tests/` are the base-R contract and regression suite,
with the golden seeded baseline. `tests/testthat/` holds the paper-fidelity
checks.
