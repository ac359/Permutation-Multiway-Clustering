# mwperm (development version)

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

# mwperm 0.2.0

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
