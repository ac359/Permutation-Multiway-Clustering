# mwperm 0.3.0

Two classes of problem are closed in this release: procedures the method paper
defines that the package never implemented, and implementation choices that made
the computed statistic differ from the paper's definition. Reported p-values and
point estimates are **unchanged everywhere**; eight seeded confidence intervals
moved, and every one of them is enumerated under *Authorized numerical changes*
below.

## New

* **`mwperm_irregular()` implements Section 6.4 (irregular designs).** The
  package previously had no implementation of it. Section 6.4 exists for two-way
  layouts where permuting *within* a cell is either invalid — the replication
  index is really time, so the errors are not exchangeable across it — or
  powerless, because `d` is constant within each cell (a dyad-level covariate).
  `mwperm_layout(L0 = )` borrowed Section 6.4's `L0` threshold but then ran the
  Section 6.3 within-cell test, which is a valid procedure but a different one,
  and left both of those cases uncovered: the package only warned "no power".

  `mwperm_irregular()` follows Section 6.4 as written. It forms the cell sizes
  and the mask `M_ij = 1{ell_ij >= L0}`, runs the biclique search
  (`find_bicliques()`, Algorithm 2) on the mask, reduces each retained cell to
  exactly `L0` observations uniformly at random (reproducibly from `seed`, with
  the usual RNG-state hygiene), and applies Procedure 2 to what remains: a row
  group on `I_q` and a column group on `J_q` per block, of common order `K + 1`,
  applied identically across the `L0` within-cell slots, so cell `(i, j)` slot
  `l` maps to cell `(pi(i), sigma(j))` slot `l`. That is the same "hold the
  third index fixed" structure `mwperm_panel()` uses for time, and it is what
  makes the test valid when the repeats are periods. Default
  `K = min over blocks of min(|I_q|, |J_q|) - 1`, capped at 199.

  Unlike `mwperm_layout()`, `d` may be constant within cells and no warning is
  issued — that is the supported case. The assumption documented is the one it
  actually needs: exchangeability across `(i, j)` within each slot, plus
  Assumption 4 on the mask; *not* within-cell exchangeability.

  Reachable as `mwperm(design = "irregular", L0 = )`, and `mwperm_check()` now
  recommends it when cells repeat and `d` is constant within every cell.

* **`aggregate = c("median", "median2")` on every front end.** Theorem 1's
  finite-sample validity is for a *single* random permutation group, so the
  p-value is exact as stated at `n_reps = 1`. The median of several dependent
  randomised p-values is a de-randomisation heuristic — the paper endorses it
  (Remark 1) and it behaves well, but it is not itself a level-alpha p-value.
  `min(1, 2 * median)` is (Rüschendorf 1982; Vovk and Wang 2020). `"median2"`
  selects it, for the reported p-value **and** for the confidence set that
  inverts it. The default is `"median"`, so no existing number changes.

* **`conf_set` field on fitted objects.** A two-column matrix of interval end
  points, one row per connected component of the confidence set. `conf_int`
  remains its hull. `ci_method` records which route produced the set:
  `"exact"`, `"grid"`, or `"bisection"`.

* **`retry_peels` argument to `find_bicliques()`** (default 4). Under the greedy
  method, a peel that returns a block below `min_block` no longer ends the
  search: the seed window slides down the degree order and retries a bounded
  number of times. A small block from the highest-degree seed rows is not proof
  that no conforming block remains elsewhere in the mask, and stopping there
  forfeited every later block too. Power only — added blocks are ordinary
  disjoint fully observed bicliques. Under `method = "exact"` a sub-floor block
  *is* proof, and the behaviour is unchanged. Verified not to change the block
  set on any shipped example; on 400 randomly generated masks with a clean block
  planted on low-degree rows it recovered extra coverage in 399 and never lost
  any.

* **A paper-to-code map** (kept with the project's development notes, not
  shipped) records every paper object —
  Assumptions 1 and 4, Procedure 1 and its three steps, Equation (10), Algorithm
  1, Procedure 2, Algorithm 2, Sections 6.1–6.4, Theorems 1 and 4 — to the
  function and file implementing it, with the reason for every deliberate
  departure.

* **`tests/golden/`** holds a seeded snapshot of every front end at default and
  non-default settings, checked by `tests/test-golden.R`. `baseline-0.2.0.rds`
  is the pre-release reference, so `Rscript tests/golden/make_baseline.R --check
  --against=baseline-0.2.0.rds` reproduces the change list below on demand.

## Corrections that change a number

* **A no-power design's p-value no longer depends on the platform's BLAS.** In
  the Section 6.3 within-cell test with `d` constant inside every cell -- the
  case `mwperm_layout()` already warns has no power -- the residualized `d` is
  constant within each cell and the permutation only reshuffles inside cells.
  The permuted statistic therefore equals the unpermuted one in exact
  arithmetic: `Dr' y_k = Dr' y` and `Dr' D_k = Dr' Dr`, so `a_k == b_k` for
  every `k` and the correct p-value is exactly 1.

  `.ipt_prepare()` recomputed `v` and `W` rather than recognising the identity.
  That sums the same terms in a different order, so they came back differing
  from `u` and `M` by about 1e-13 -- and that rounding noise, not the data, then
  decided the `a_j <= b_k` comparisons. The reported p-value was whatever the
  local BLAS happened to produce: 1 on this project's macOS build, 0.75 or 0.5
  on the Linux and Windows CI runners, from bit-identical inputs.

  When the gather leaves the residualized regressor unchanged
  (`identical(Dr[g, ], Dr)`), the package now asserts the identity instead of
  recomputing it. Degenerate fits report p = 1 on every platform; every
  non-degenerate fit takes the same path as before and is bit-identical, which
  the golden baseline confirms. Only designs the package already warns are
  powerless are affected -- no valid inference changes -- but a seeded p-value
  did move on non-macOS platforms, so it is recorded here.

* **The confidence set is now computed exactly, not by bisection.** Procedure 1,
  step 3 defines `CI = {b : pval(b) > alpha}`. `.invert_ci()` approximated it: it
  assumed a single interval, bracketed outward, bisected to a tolerance, and
  patched disconnected cases with a hull. For a single coefficient the set is
  available in closed form from the cached cross products. The statistics
  `a_j(b) = |u_j - M_j b|` and `b_k(b) = |v_k - W_k b|` are piecewise linear, so
  the p-value is a step function of `b` whose only jumps solve
  `|v_k - W_k b| = |u_j - M_j b|`, at `b = (u_j - v_k)/(M_j - W_k)` and
  `b = (u_j + v_k)/(M_j + W_k)`. The package now pools those roots across
  repetitions and evaluates the aggregated p-value at every root and inside every
  interval between consecutive roots, plus one point beyond each end. No
  bracketing assumption, no tolerance, and a disconnected set is reported as its
  components rather than replaced by its hull. Above a candidate budget
  (about `2 K^2 * n_reps`, default cap 200,000, settable with
  `options(mwperm.ci_exact_budget = )`) it falls back to the old
  bracket-and-bisect path and says so in `$note`.

* **One rep-aggregation rule now drives every path.** Three code paths
  aggregated `n_reps` differently and two of them were not inversions of the
  same function as the reported p-value: the reported p-value took the median
  across reps, but the default interval took the **median of the per-rep end
  points**, and the explicit-`grid` path had (before 0.2.0) accepted a point if
  *any* rep accepted it. The confidence set is now defined once, as
  `{b : median_r pval_r(b) > alpha}`, and the exact path, the grid path, the
  bracketing fallback and the joint region all compute that one set. The
  reported p-value goes through the same function. A value inside the reported
  interval can therefore no longer be rejected by the test, and vice versa.

  This is what moves the eight intervals below. Median-of-end-points was
  systematically **narrower** than the inversion of the median p-value, so the
  old default intervals were slightly too narrow; every changed interval is
  wider.

* **Sub-seed derivation is overflow-safe and collision-free.** `.sub_seed()`
  computed `rep_seed * 1000L + j`. Two problems: an integer `seed` above roughly
  2.1e6 overflowed to `NA` and `set.seed(NA)` failed; and with 1000 or more cells
  (layouts) or 250 or more blocks (missing and irregular designs) two `(rep, j)`
  pairs collided and shared a relabelling. A collision never broke validity —
  each rep's test is exact regardless, the seed only picks the random
  relabelling — but the reps were not independent, so the median was averaging
  fewer effective draws than it appeared to. Fixed in two value-preserving
  steps: the arithmetic is done in double, and the stride is widened **only** by
  callers whose `j`-range would collide. Every design that was already
  collision-free keeps exactly the seeds, and the results, it had.

## Corrections that change no number

* **Confidence sets are 2.6-3.4x faster, with every result bit-identical.**
  The exact set introduced in this release evaluates the aggregated p-value at
  every breakpoint and every cell between breakpoints, so for a 40x40 dyadic
  fit at defaults it aggregates a 60,841 x 10 matrix of per-repetition
  p-values. Two hot spots in `R/engine.R` were rewritten:
  `.agg_pvals()` now uses a vectorised `.row_median()` instead of
  `apply(P, 1L, stats::median)` -- which was making one R-level `median()` call
  per row, profiled at ~49% of a whole fit -- and `.pval_matrix()` streams its
  inner block instead of materialising two `K x chunk` matrices, a logical
  matrix and a recycled vector.

  Nothing about the method changed: the statistic, the minorisation, the
  aggregation rule, the breakpoint enumeration and the set definition
  `{b : agg_r p_r(b) > alpha}` are all untouched. Bit-identity was the
  precondition, because the aggregated value is compared to `alpha` with `>`
  and per-repetition p-values sit exactly on the grid `j/(K+1)`, so ties at
  `alpha` are common. Note for anyone revisiting this: `(lo + hi)/2` is *not*
  bit-identical to `mean(c(lo, hi))` -- `mean()` applies a second-pass LDOUBLE
  correction, and LDOUBLE is 80-bit on x86 but 64-bit on arm64 -- so
  `.row_median()` obtains order statistics by selection and calls `mean()`
  itself, once per distinct pair.

  Verified: `tests/golden/make_baseline.R --check` reproduces, all 19 test
  files pass, and a 22-fit battery spanning all six designs, `n_reps` 1 to 15,
  both `aggregate` rules, explicit `grid`, `d > 1`, non-zero `beta_null` and
  `n_cores = 4` is bit-identical fit for fit. Measured: `trade_dyadic` at
  defaults 1.70s -> 0.52s, `trade_panel` 0.61s -> 0.32s. Fits without a
  confidence interval are unaffected.


* **The rejection floor is now aggregation-aware, so `aggregate = "median2"`
  reports its resolution limit instead of an unbounded interval.** `"median2"`
  reports `min(1, 2 * median)`, so its smallest attainable p-value is `2/(K+1)`,
  not `1/(K+1)`. The engine gated both the "cannot reject at this level" note
  and the confidence set on `1/(K+1) > alpha` regardless of the rule. At
  `K = 19` and `alpha = 0.05` that test saw `0.05 <= 0.05` and stayed silent
  while `"median2"` could never return below `0.10`: on a 20x20 design with an
  overwhelming effect it returned `p = 0.10` and `conf_int = [-Inf, Inf]` with
  no note, where the default rule rejected and returned a finite interval. The
  engine now computes the floor from the aggregation rule and gates on that, and
  the note names it (`2/(K+1)`) and asks for `K + 1 >= 2/alpha` — 40 levels in
  the smallest permuted dimension at `alpha = 0.05`, twice what `"median"`
  needs. The floor is also stored as the new `p_floor` field and is what
  `confint()` cites when no set was computed; `$resolution` is unchanged and
  remains the p-value grid step `1/(K+1)`, as printed.

  **Effect on output.** Under the default `aggregate = "median"` nothing moves:
  the floor equals the grid step, so every gate, note and interval is what it
  was, and the golden baseline reproduces. Under `"median2"` on a design too
  coarse for the level, a note now fires and the unbounded interval becomes
  `NULL` with that explanation — the convention the guard already used
  everywhere else. No finite interval changes under either rule.

* **The confidence set is documented as the *closure* of `{b : pval(b) > alpha}`,
  which is what it has always been.** The exact route builds each component from
  a maximal run of accepted atoms of the step function and reports the bounding
  breakpoints, so a component that begins or ends strictly inside an open cell
  is reported with the breakpoint just outside it — a value the test rejects.
  On the `trade_dyadic` anchor fit, `conf_int[1] = -1.251351` has `p = 0.0500`
  (rejected at `alpha = 0.05`) while `p(-1.251351 + 1e-12) = 0.0625`. The
  convention is the usual one for a discrete p-value and errs outward: no value
  the test accepts falls outside the reported set, and every point strictly
  inside a component is accepted. `README.md` and `?confint.mwperm` previously
  claimed an exact if-and-only-if, which does not hold at the boundary; both now
  state the closure relation, and `tests/lower-level-tests/test-exact-ci.R`
  probes the breakpoints
  themselves and pins it. Numbers are unchanged — only the claim was wrong.
  (The `"grid"` and `"bisection"` routes report attained, accepted end points.)

* **Every permutation is asserted to be a bijection.** `.build_obs_perms()` and
  the block-diagonal builder in `mwperm_missing()` map permuted cells back to
  observation indices with `match()` on a mixed-radix cell code. `match()`
  returns the first hit, so two observations sharing a cell code would make the
  gather vector many-to-one: the statistic would be computed on duplicated rows,
  with no NA, no warning, and no way to notice. Both builders now reject
  duplicate cells at entry, and every returned gather vector is checked to be a
  permutation of `seq_len(N)` before it is used. The checks are O(N) against
  O(N p^2) of linear algebra, so they run unconditionally. This converts a
  silent wrong answer into an error; it cannot change a correct result.

* **A permutation that annihilates `d` is now an error.** If the residualized
  `d` is numerically zero for some permutation `k`, that permutation's statistic
  is rounding noise rather than data, and the a/b comparison is decided by float
  error. It was silently zeroed. It now stops, naming the design. The
  pre-existing case where `d` is constant or collinear with `x` — where *every*
  slice is degenerate and `p = 1` is the exact answer — still warns and returns
  `p = 1` as before.

* **The projection dimension claim is corrected.** The paper writes
  `V_k` in `R^{N x (N - 2p)}`; the code projects onto the orthogonal complement
  of the column space of `[X | X_k]`, whose dimension is
  `N - rank([X | X_k])` and is strictly larger whenever a nuisance column is
  permutation-invariant — always, through the intercept, and substantially in
  `mwperm_panel(time_fe = TRUE)`, where the period dummies are permuted among
  themselves. The code's choice is the correct one and `N - 2p` is not well
  defined under rank deficiency, so the claim was fixed, not the code, in
  `README.md`, `?mwperm_dyadic`, `?mwperm_panel` and the corresponding `.Rd`
  files.

* **The Section 6.3 / 6.4 attribution of `L0` is corrected.** `L0` stays in
  `mwperm_layout()` — dropping thin cells raises the attainable `K` and so the
  p-value resolution, which is a legitimate power lever — but the roxygen block
  and `man/mwperm_layout.Rd` now state that the threshold comes from Section
  6.4, that `mwperm_layout()` applies it to the Section 6.3 within-cell test,
  and that `mwperm_irregular()` is the Section 6.4 procedure.

* **`?confint.mwperm` documents the set.** A new section states the definition
  the package inverts, the aggregation rule across repetitions, that the set is
  computed exactly for a single coefficient, and what `conf_set` and `ci_method`
  contain.

* **`?mwperm_dyadic` notes what the default `K` does.** `K = min(n_row, n_col) - 1`
  maximises p-value resolution and makes `K + 1` equal the cluster count, so
  Algorithm 1 produces a single block spanning all clusters. Theorem 2's power
  result is stated for a *fixed* `K` with the cluster count divisible by
  `K + 1` — a different regime. Validity is unaffected either way, so the
  default stands; the note says so rather than leaving the reader to assume the
  power theorem describes the default.

* **`README.md`** documents `mwperm_irregular()` in the design-choice table and
  the function map, describes the exact confidence set and the single
  aggregation rule, and records that the test suite is tracked and ships in the
  tarball (`tests/` is *not* excluded from the public branch, contrary to an
  older note).

* **Citations verified against the paper's title page.** The method paper is by
  Wenxuan Guo, Panos Toulis and Yuhao Wang, arXiv:2601.08610. Every `.Rd` file,
  roxygen block, `README.md`, `NEWS.md`, the vignette and `inst/CITATION`
  already agreed on this; `DESCRIPTION` gave no initial for Guo and now reads
  "Guo, W., Toulis, P. and Wang, Y. (2026)". (F. R. Guo, whom the method paper
  separately cites as Guo and Shah (2025), is a different statistician and does
  not appear anywhere in the package.)

* **`DESCRIPTION`.** Version 0.3.0. `grDevices` was already declared in
  `Imports` and imported in `NAMESPACE`; verified against the source (`R/plot.R`
  uses `chull`, `adjustcolor` and the device openers) and left as is.

* **`README.md`: two shown outputs were still the 0.2.0 numbers**, missed when
  the exact confidence set changed eight seeded end points earlier in this
  release. The quick-start block showed `[-1.246, -0.5486]` where the package
  prints `[-1.246, -0.5485]`, and the panel example under *Extensions* showed
  `[0.4441, 0.8776]` where it prints `[0.442, 0.8803]`. Both are corrected; the
  package's own output never changed, only the transcript of it. Every shown
  output in `README.md` was re-run against the installed package.

* **"Choosing the right design" now gives the formal model for each row.** The
  table gained a column with the error decomposition each design is exact
  under -- $\varepsilon_{ij} = \eta_i + \xi_j + u_{ij}$ for the two-way case,
  $\varepsilon_{ijt} = \eta_i + \xi_j + \zeta_t + u_{ijt}$ with $\zeta_t$
  arbitrary for the panel, and so on -- plus a column saying which labels are
  permuted. A new subsection states the invariance each design actually needs
  as a display equation (InvA, InvB, within-cell, blockwise InvB, and
  Assumption 4 on the mask), and a four-question checklist walks a reader who
  is still unsure to a function. The invariance statements that were duplicated
  under *Extensions* were removed from there, so each condition is now stated
  in exactly one place.

* **The test suite is reorganised around the public API, and documented.**
  `tests/` was organised by concern -- `test-edgecases.R`, `test-paths.R`,
  `test-fixes.R` -- so there was no file to open to find out how a given front
  end was tested, and three files accumulated most of the assertions for all
  five designs. The top level now holds **one file per user-facing entry
  point** (`test-dyadic.R`, `test-threeway.R`, `test-panel.R`,
  `test-layout.R`, `test-missing.R`, `test-irregular.R`, `test-formula.R`,
  `test-main.R`), plus the cross-cutting contracts (`test-validation.R`,
  `test-methods.R`, `test-equivariance.R`) and the seeded gate
  (`test-golden.R`). Tests of the machinery underneath -- the permutation
  group, the gather vectors, the projector, the p-value, the aggregation rule,
  the confidence-set routines -- moved to `tests/lower-level-tests/`, driven by
  `test-lower-level.R` so `R CMD check` still runs every one of them.
  Assertions shared by several files moved to `tests/helpers/assertions.R`, and
  `tests/README.md` documents the layout and where a new test belongs. **No
  assertion was dropped**; the suite is the same coverage, reorganised, plus
  new per-design checks that had no home before -- the default `K` each design
  derives, dense recoding of arbitrary cluster labels, the panel test's
  invariance to an arbitrary common time trend, and row-order invariance for
  every design.

## Authorized numerical changes, in full

Reported **p-values and point estimates are identical everywhere**. Biclique
block sets on the shipped examples are identical. The complete list of moved
numbers, reproducible with
`Rscript tests/golden/make_baseline.R --check --against=baseline-0.2.0.rds`:

| Design and argument regime | Field | 0.2.0 | 0.3.0 | Why the old number was wrong |
|---|---|---|---|---|
| `mwperm_dyadic()` on `trade_dyadic`, defaults (`n_reps = 10`, `seed = 1`); identically via `mwperm()` and `mwperm_formula()` | `conf_int` | `[-1.244156, -0.542957]` | `[-1.251351, -0.534468]` | Median of per-rep end points, not the inversion of the median p-value. Too narrow by 2.2% of width. |
| `mwperm_panel()` on `trade_panel`, defaults (`n_reps = 10`, `seed = 1`); identically via `mwperm()` | `conf_int` | `[0.444144, 0.877645]` | `[0.442022, 0.880272]` | Same cause. Too narrow by 1.1% of width. |
| `mwperm()` README quick start on `trade_dyadic` (`n_reps = 15`, `seed = 1`) | `conf_int` | `[-1.246414, -0.548561]` | `[-1.246471, -0.548496]` | Same cause plus the removed bisection tolerance; 0.017% of width. |
| `mwperm_dyadic()` on `trade_dyadic`, `beta_null = -1`, `n_reps = 9`, `seed = 1` | `conf_int` | `[-1.236987, -0.534473]` | `[-1.237042, -0.534468]` | Bisection tolerance removed; 0.008% of width. |
| `mwperm_panel()` on `trade_panel`, `beta_null = 0.5`, `n_reps = 9`, `seed = 1` | `conf_int` | `[0.442033, 0.880264]` | `[0.442022, 0.880272]` | Bisection tolerance removed; 0.004% of width. |

Every changed interval is **wider**. Not represented in the table, because no
shipped example triggers them, but authorized and possible:

* **Any single-coefficient `conf_int` with `n_reps > 1`** can move, for the
  reason above. The direction is not guaranteed in general, but median of end
  points is narrower than the inverted median whenever the per-rep sets are
  nested or nearly so, which is the usual case.
* **Any `conf_int` produced with an explicit `grid`** now reports the components
  in `conf_set`; the interval itself is unchanged from 0.2.0 (the grid path
  already inverted the median).
* **A disconnected confidence set** now reports its components in `conf_set`;
  `conf_int` is still the hull, so the reported interval is unchanged, but the
  accompanying note is more specific.
* **Layouts with 1000 or more occupied cells, and missing or irregular designs
  with 250 or more blocks**, get different seeded output, because the sub-seed
  stride is widened for exactly those designs. Below those sizes nothing
  changes. No shipped example or documented number is in that regime.
* **`find_bicliques(method = "greedy")` on masks where a peel returns a
  sub-floor block** may now return more blocks, which changes the retained cells
  and hence every downstream number for such a fit. `retry_peels = 0L` restores
  the old behaviour. No shipped example is affected.

# mwperm 0.2.0

* **Documentation rewritten against verified behaviour.** `README.md` was
  reorganised so that a reader reaches a correct first call before meeting any
  mathematics, and every code block and shown result in it is now real pasted
  output rather than an abbreviation. Several claims were corrected: the
  package was described as "valid with as few as ~20 clusters per dimension",
  which conflated validity with resolution and understated the guarantee (the
  test is exact at *any* number of clusters; 20 is what is needed to be able
  to reject at the 5% level); `mwperm_formula()`, `coef()` and `nobs()` were
  missing from the function list; `grDevices` was missing from the stated
  imports. New sections cover choosing a design, reading the output,
  reproducibility, performance, and limitations.

* **`print()` now states the p-value resolution.** The p-value is exact but
  discrete -- it can only take multiples of `1/(K+1)` -- and nothing in the
  output said so. A fit reporting `p-value = 0.025` at `K = 39` was reporting
  the *smallest value that fit could ever produce*, which is easy to read as
  a precise number; conversely a design with `1/(K+1) > alpha` cannot reject
  at `alpha` however large the effect. Every fit now prints, for example,
  `Resolution   : p-values are multiples of 1/40 = 0.025`.

* **Clearer output for joint (multi-coefficient) tests.** The `H0` line marks
  the null as joint over all coefficients (a scalar `beta_null` recycled to
  several printed as a bare `0`, reading as a single-coefficient test), and
  the footer now says the printed brackets are one confidence region's
  marginal extent rather than separate per-coefficient intervals.

* **Better guidance in notes, warnings and errors.** The coarse-resolution
  notes advised "increase K", which is not actionable because K is capped by
  the design; they now give the concrete requirement (at least `1/alpha`
  levels in the smallest permuted dimension) and state that the p-value
  itself remains exact. `confint()` no longer suggests refitting with
  `conf_int = TRUE` when `conf_int` was already `TRUE` and the resolution was
  the binding constraint. The exact-biclique budget warning now says validity
  is unaffected (a sub-maximal block costs power only). `mwperm_missing()`
  reports discarded cells as the mechanism of Procedure 2 rather than as a
  bare percentage. Layout fits label their two counts as cells and replicates
  instead of "clusters", and `K = NULL` errors name the admissible range.

  These are wording, formatting and documentation changes only: no p-value,
  estimate, confidence limit, `K`, or note-firing condition changes. Verified
  against a 29-fit regression battery in which 19 objects are `identical()`
  outright and the other 10 differ only in the free-text `note` field.

* **Speed: the permutation builder is 2.2-2.6x faster, with bit-identical
  results.** `.build_obs_perms()` recomputed the mixed-radix bases on every
  one of the K+1 permutations (through `apply(coords, 2L, max)`, which copies
  the whole matrix), and materialized a permuted coordinate matrix each time.
  The bases cannot change under a permutation -- each image vector is a
  bijection of that dimension's ids -- so they are now computed once and
  reused, and the cell coding is inlined over the permuted coordinates
  instead of building the matrix. End-to-end this is 1.17-1.24x on the
  dyadic, panel and three-way designs (panel with 200 clusters and 10
  repetitions: 110 s to 93 s); the layout and missing-data designs use
  different builders and are unchanged.

  **No result changes.** `max` is a comparison reduction, so no arithmetic
  was reassociated, and the inlined coding uses the same bases, place values
  and accumulation order as before. Verified against the previous
  implementation pasted in verbatim: cell codes identical on 6 of 6 coordinate
  layouts, gather vectors identical on 8 of 8 designs (including the sparse
  `match()` branch), and a 29-fit regression battery spanning all five
  designs identical under `identical()`.

* **Packaging fix.** `DESCRIPTION` gave the maintainer (`cre`) role to all
  three authors; R permits exactly one, so the package failed to install --
  `R CMD INSTALL`, `build`, `check` and every CI matrix leg stopped at the
  DESCRIPTION parse step before any R code ran. The `cre` role is now held by
  the maintainer of record alone; all three authors keep `aut` and `cph`,
  matching `LICENSE`. No behaviour or result changes.

* **Bug fix: `d` supplied as a value rather than an expression.** Passing the
  covariate of interest through `do.call()` -- a normal way to drive the
  package programmatically -- made every design-specific front end fail with
  an opaque `'names' attribute [N] must be the same length as the vector [1]`,
  because `deparse(substitute(d))` returns one element for a symbol but many
  for a value. `mwperm()` did not error but labelled the coefficient with the
  entire deparsed data vector, which then appeared in `print()`,
  `summary()$term`, `confint()` row names and plot axes. Both paths now fall
  back to the generic label `"d"` when the deparse is not a single element.
  Calls that already worked are unaffected: their labels, p-values and
  confidence limits are unchanged.

* **Documentation corrections.**
  - `?mwperm_layout` attributed the `L0` balancing threshold to Section 6.3 of
    Guo, Toulis and Wang (2026); `L0` is defined in Section **6.4**. The page
    now cites 6.4 and states plainly which part of that procedure is
    implemented: `L0` performs 6.4's balancing step and then applies 6.3's
    within-cell test, so a covariate that is constant within cells still
    yields a powerless test and replicates that are really time periods are
    still not permutable.
  - `?mwperm`'s `grid` argument still described the pre-0.2.0 behaviour (the
    interval as the range of grid points not rejected). It now matches the
    shipped median-aggregated semantics already documented on the
    design-specific pages.
  - The package page's author block now records the copyright holders.

* **Performance (no change to any result).** The per-permutation kernel was
  rebuilt three times. It now residualizes only `d` on the stacked design
  `[X | X_pi]` -- `V_k V_k'` is symmetric and idempotent, so the outcome and
  the permuted columns enter the statistic as plain inner products against
  the residualized `d`, and the three extra columns the previous code pushed
  through the same back-substitution produced no new information. It then
  drops the per-permutation QR entirely: because `X_pi` is a row permutation
  of `X`, the two diagonal blocks of the stacked design's Gram matrix are
  both `X'X`, so `X'X` is formed once per fit and only the cross block
  `X'X_pi` -- one `crossprod` -- depends on the permutation. The
  block-permutation builder behind `mwperm_missing()` also gained the
  position-table cell lookup already used elsewhere. Finally, that cross
  block is now spelled `t(X) %*% X_pi` rather than `crossprod(X, X_pi)`:
  identical sums over the same reduction index, but reaching a BLAS kernel
  that accumulates each entry with independent accumulators instead of one
  serial dependency chain. The block-permutation builder additionally
  composes each block's maps at block-label rather than cell length. A
  representative gravity fit (`mwperm_missing`, N = 9,248, K = 67,
  `n_reps = 15`) drops from 11.7 s to 2.8 s serial. P-values, per-rep
  p-values, estimates and intervals are unchanged -- verified `identical()`
  across all five designs, on all 28 interval endpoints of the gravity
  application, and against the authors' own implementation on a fixed
  permutation group. The relative gain of the last step depends on the BLAS
  in use.

* **Bug fix (changes reported intervals).** With an explicit `grid` and
  `n_reps > 1`, the single-coefficient confidence interval retained any value
  accepted in *any* repetition (a union, i.e. inversion of the maximum
  p-value across reps). It now inverts the **median** p-value across reps,
  matching the reported p-value, the joint-region path, and Remark 1 of Guo,
  Toulis & Wang (2026). Grid-mode intervals are now narrower and no longer
  grow with `n_reps`. Intervals whose acceptance region reaches the edge of
  the supplied `grid` are now reported as unbounded on that side rather than
  silently truncated.

Audit-driven release: the changes below implement the fix plan from the
2026-07 verification audit (finite-sample validity re-verified throughout;
no change to the test statistic, the permutation group, the minorized
p-value, or any seeded result unless explicitly noted).

* **License:** the package is now released under the MIT license
  (`LICENSE`); previous versions carried a placeholder `Proprietary` tag.
* **Default `n_reps` is now 10** (was 1) in every test function and in
  `mwperm()`. A single run's p-value depends on the random relabelling --
  the audit documented a real seed lottery -- while the package's own
  documentation has always recommended median aggregation; ten repetitions
  cost fractions of a second on typical designs and were measured
  conservative (never anti-conservative) at every level. **This changes
  default-argument seeded results** (calls with an explicit `n_reps` are
  unaffected): reported p-values and confidence limits are now medians over
  10 runs. Set `n_reps = 1` to reproduce pre-0.2.0 defaults.
* `DESCRIPTION` gains `URL`, `BugReports`, `Language`, `Date`; a citation
  entry ships as `inst/CITATION` (`citation("mwperm")`).
* **New formula interface:** `mwperm_formula(y ~ d | x, data, index, ...)`
  (the nuisance part is optional). Both formula parts are standard formula
  algebra via `model.matrix()`, so transformed terms and factors work; the
  result is identical to the data interface with the same seed (pinned by
  tests). New accessors `coef()` (the OLS estimate, named) and `nobs()`.
* **New `permute` argument on `mwperm_missing()`** (and passed through by
  `mwperm()`): `"rows"` or `"cols"` permutes only that dimension. One-sided
  permutation applies a subgroup of the invariance group, so the test stays
  exactly valid under the same exchangeability assumption, and `K + 1` is
  then capped by the *permuted* block side alone -- blocks need `min_block`
  clusters on the permuted side but as few as one on the other, so designs
  whose fully observed blocks are short in one dimension (down to a single
  fully observed column) can reach resolutions the two-sided test cannot,
  at some cost in power. `find_bicliques()` correspondingly accepts a
  length-2 `min_block = c(rows, cols)` and, when the area-maximal block
  violates such an asymmetric floor, retries the greedy growth under the
  floor (a tall thin block never maximises area on a dense mask, so without
  the retry the tall blocks the floor asks for would not be found). A scalar
  `min_block` keeps the historical both-sides floor bit-identically, and the
  default `permute = "both"` path is unchanged. Size control and power of
  the one-sided test were verified by simulation (2,000-run null cells,
  including the constrained-retry path; no super-uniformity violation).
* `R CMD check --as-cran` passes with no errors and no package-level
  warnings. Two fixes were needed: the `coef()`/`nobs()` S3 methods (new in
  this release) registered without importing their generics from `stats`, so
  loading the namespace in isolation failed with `object 'nobs' not found`
  (ordinary `library(mwperm)` masked it because `stats` is always attached) --
  both generics are now imported; and `n_cores` is now additionally capped by
  `getOption("mc.cores")` when set, so the standard core-throttling option
  (and R CMD check's core limit) is honoured. The `.github/` CI directory is
  excluded from the build tarball.
* **Replication material** ships in `inst/replication/`: numbered,
  self-contained Monte-Carlo scripts (base R + `mwperm` only) that reproduce
  the package's headline claims -- finite-sample size control versus the
  invalid classical OLS *t*-test, power rising to 1, confidence-interval
  coverage, the trending-panel negative control behind the panel-by-default
  policy, and the one-sided `permute` option -- with a `make.R` driver,
  cached/resumable runners, reference outputs under `expected/`, and a
  recorded `sessionInfo()`. See `system.file("replication", package =
  "mwperm")`.
* **Style and documentation-completeness pass** (no behaviour change; the
  full suite and every seeded reference value verified bit-identical): all
  code lines are now <= 80 characters with no compound/trailing semicolons;
  source and tests are pure ASCII; every exported function's help page has
  runnable examples, references, and cross-references, and the package
  spell-checks clean against a shipped `inst/WORDLIST`.
* Shipped-test coverage of `R/` is 94% (audit gate: >= 90%; it was 81% at
  audit time): new `tests/test-paths.R` exercises the design-diagnosis
  printer, forced-design validation, name resolution, `L0` balancing, the
  explicit-grid and joint-region confidence paths, a constructed
  disconnected acceptance set (the island guard's hull), and every figure
  export device. GitHub Actions CI runs `R CMD check --as-cran` on
  ubuntu/macOS/windows across devel/release/oldrel-1 plus a coverage job.
* The test suite (9 base-R files, ~2 s under `R CMD check`) is now tracked
  in the repository and ships with the package; internal development
  material is excluded from builds via `.Rbuildignore`.
* **Detection safety** (`mwperm()`/`mwperm_check()`): the time role is never
  assigned anti-conservatively in silence. A name-based time assignment that
  the values do not corroborate now warns (a column merely *named* "period"/
  "year" may be a cluster; permuting the true time dimension can over-reject
  badly -- the audit measured 87% rejection at a nominal 12.5% on such a
  case). Forcing `design = "threeway"` when an index looks time-like also
  warns, and the ambiguous-fork message now states correctly that the panel
  default protects only if the *permuted pair* is exchangeable. All
  detection warnings are raised as R warnings at fit time (previously they
  were only recorded on the returned object's `note`). Detection *choices*
  are unchanged -- same data, same seed, same result.
* A factor `y` is now rejected with an informative error in every front end;
  it was previously coerced silently to its internal level codes (`d` and
  `x` were already protected).
* **Degenerate `d` is now deterministic:** when `d` has no variation after
  partialling out the nuisance covariates (constant, or collinear with a
  column of `x`), the test warns that beta is unidentified and returns
  p = 1 exactly -- the exact-arithmetic answer, since the residualized
  statistic is identically zero and the minorized p-value is 1. Previous
  versions let ~1e-16 rounding noise decide the permutation comparisons,
  giving an arbitrary, BLAS-dependent p-value with no warning. Only
  degenerate fits are affected; all seeded reference values are unchanged.
* **Validation hardening:** `K` is validated as a single integer >= 1 with
  an error naming the argument (previously `K = NA`/vector `K` surfaced raw
  R errors, a fractional `K` truncated silently, and `K = 0` blamed the
  data); `seed` must be `NULL` or a single finite number, and a seed too
  large for the rep/sub-seed scheme (|seed| above ~2.1 million) now errors
  by name instead of via `set.seed(NA)`'s cryptic message -- in-range seeded
  results are bit-identical; an `NA` in the layout `rep` identifier is
  rejected (it was silently ranked last within its cell);
  `mwperm_missing()` validates the full supplied data before the biclique
  step discards cells; `confint(parm = )` now subsets by coefficient name
  or position (it was accepted and ignored).
* `mwperm_missing()` now says so in its `note` when the smallest selected
  block caps `K` below the level's resolution (no rejection attainable at
  `alpha`), naming the block as the cause and `min_block` as the remedy;
  `?mwperm_missing` documents that the smallest block is the binding
  constraint. The block-diagonal permutation builder was extracted to a
  named internal so the shipped group-closure tests exercise the exact
  production code (no behaviour change).
* **Documentation is now single-source:** `man/` and `NAMESPACE` are
  generated by roxygen2 from the comments in `R/` (they were hand-written
  and had drifted). All hand-written surplus was ported into the roxygen
  comments first -- examples for `mwperm_layout()`/`mwperm_threeway()`/
  `plot.mwperm()`/`mwperm_save()`, the full plot style-element list,
  `\seealso` sections, dataset pages (now with `\source`) -- so no rendered
  content was lost; the front-end pages gain the full inherited argument
  documentation their hand versions abbreviated.
* **Parallel path fixes:** on non-fork platforms (Windows) the engine now
  creates ONE worker cluster per fit and reuses it across repetitions --
  previously a fresh PSOCK cluster was spawned inside every repetition,
  which made `n_cores > 1` about 3x *slower* than serial under the default
  settings; `n_cores` beyond the detected core count is now clamped with a
  warning (it was silently oversubscribing). Parallel results remain
  bit-identical to serial.
* **Performance** (bit-identical output; verified against the previous
  implementation on seeded fits across all five designs, plus the full test
  suite and all seeded reference values): layout fits are ~10x faster (the
  per-cell bookkeeping in the permutation builder and the within-cell index
  construction are now vectorized -- 5.3 s -> 0.44 s on a 5,000-cell
  layout); every fit gains ~20% from a fused single-call residualization
  (`.lm.fit`, same pivoted-QR family as `qr.resid`, so rank-deficient
  stacked designs are handled identically); complete-array gather vectors
  are built by an O(N) position-table translation instead of per-element
  hashing.
* **Documentation honesty sweep** (audit findings): the exact status of
  median aggregation over `n_reps` is stated (each repetition is valid on
  its own; no finite-sample theorem covers the median; measured uniformly
  conservative; the 2x-median rule is the provable fallback);
  reproducibility of seeded results is documented as conditional on
  `RNGkind()` and, for character ids, the collation locale (validity is
  unaffected); near-collinear nuisance columns are documented as silently
  dropped at the QR tolerance; `?mwperm_layout` no longer presents a
  shared-across-cells replicate effect as covered by the within-cell
  invariance argument; `?find_bicliques` states that `"exact"` maximises
  each block in turn (not total coverage) and that `min_block` is floored
  at 2; the panel `N > 2p` feasibility requirement is documented as
  conservative under `time_fe = TRUE`; `?build_perm_set` records that group
  closure is certified by the algebraic test suite, not by simulation;
  parallelism docs now say only the seeded repetition axis pays.

# mwperm 0.1.1

Publication-standard plotting and provenance-labelled return values (no
change to any test, seed, or p-value).

* **Breaking:** the `summary()` data frame's columns now name the source of
  each quantity, matching the printed "OLS estimate" / "IPT CI" labels:
  `estimate` -> `ols_estimate`, `se_naive` -> `ols_se_naive`,
  `conf_low` -> `ipt_ci_low`, `conf_high` -> `ipt_ci_high` (`term` and
  `p_value` are unchanged; `p_value` is the IPT permutation p-value). No
  numeric value changed. Code that indexes the old column names must be
  updated; no deprecated aliases are provided.
* `confint()` still returns the percentile-labelled matrix the generic
  promises (`"2.5 %"`/`"97.5 %"`), but now records its provenance in a
  `"method"` attribute, `"IPT (inverted permutation test)"`. The `\value`
  documentation of every test function now states which fields of the
  returned object are the OLS estimate/naive SE (`estimate`, `se_naive` --
  the SE is only the centre/scale of the confidence-set search) and which
  are the IPT confidence set (`conf_int`, `conf_region`, `conf_box`); the
  field names themselves are unchanged.
* `plot.mwperm()` gains a `type` argument: `"coef"` (OLS estimate against the
  inverted-test confidence set -- an interval with end caps for one
  coefficient, a forest of the joint region's marginal extents for several),
  `"region"` (joint confidence region, two coefficients), `"stability"` (the
  Monte-Carlo p-value diagnostic, restyled), and `"all"`. **The default
  figure changed:** `type = "auto"` now draws the flagship `"coef"` figure
  whenever a confidence set is stored, falling back to `"stability"`
  otherwise (previously a p-value histogram/barplot was always drawn).
* All figures share one style layer: Okabe-Ito colourblind-safe palette with
  marks/line types differing as well as colour (grayscale-legible), no
  top/right spines, and a standard annotation block (design, cluster counts,
  N, resolution 1/(K+1), the null, the p-value, and the decision at alpha).
  Style elements can be overridden by name, e.g.
  `plot(fit, col_estimate = "black")`; unknown `...` arguments now warn and
  are ignored instead of crashing the underlying graphics calls.
* `"null"` and `"profile"` types are reserved (permutation-null and
  test-inversion p(b) figures); they need fit-time storage that this version
  does not retain and currently fall back to the default figure with a
  message.
* New `mwperm_save()` writes any figure at journal dimensions (single-column
  3.5 in / double-column 7 in, >= 300 dpi; pdf/png/tiff/jpeg).
* `grDevices` added to Imports (still base R only). Plotting never mutates
  global graphics state (`par()` is restored on exit).

# mwperm 0.1.0

Initial release.

* Finite-sample-valid invariant permutation tests (Guo, Toulis & Wang, 2026)
  for regression coefficients under multi-way clustering:
  `mwperm_dyadic()`, `mwperm_threeway()`, `mwperm_panel()` (arbitrary common
  time trend), `mwperm_layout()` (replicated two-way layouts, with `L0`
  balancing), and `mwperm_missing()` (incomplete arrays via fully observed
  bicliques, greedy or exact solver).
* Unified entry points: `mwperm()` auto-detects the design from the
  clustering structure and dispatches (with identical results to the direct
  call); `mwperm_check()` prints the diagnosis without running anything.
  Assumption-dependent forks (panel vs three-way; layout vs suppressed
  panel) default to the choice valid under the widest set of error processes
  and are announced with override instructions.
* Confidence sets by test inversion: an interval for one coefficient
  (with a guard that detects disconnected acceptance regions and widens to
  their hull with a note), a grid-based joint region for several.
* Opt-in parallelism via `n_cores` on all test functions (forked workers on
  Unix, PSOCK on Windows); results are identical to serial runs.
* S3 methods `print()`, `summary()`, `confint()`, `plot()`; synthetic
  example data `trade_dyadic` and `trade_panel` (generated reproducibly by
  `data-raw/make_data.R`).
