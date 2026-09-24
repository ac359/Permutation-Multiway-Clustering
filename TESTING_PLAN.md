# TESTING_PLAN.md: paper-fidelity test suite and source-comment pass

Working checklist for the session that adds (a) a testthat suite showing that
`mwperm` computes what Guo, Toulis & Wang (2026) and Wen, Wang & Wang (2025)
define, and (b) source comments that map every step of the method to its
code. It is kept current so the work can resume after a context reset. The
file is listed in `.Rbuildignore`, so it is not in the tarball.

**Sources.**

* **IPT paper:** Guo, Toulis & Wang (2026), arXiv:2601.08610v1.
  It is the authority on the mathematics. The sign-flip assumption is the
  revised paper's Assumption 2 ("double sign symmetry"). That text is not in
  v1; the maintainer supplied it as an excerpt on 2026-09-24.
* **RPT paper:** Wen, Wang & Wang (2025), *Ann. Statist.* 53(2).
  Its PDF text layer is letter-shifted, so it was read from the page images.
* **JSS draft:** the 2026-09-18 build, which describes package 0.4.1. It is
  the authority on the intended API and on user-facing behavior.
  **Section numbering:** in this build, §6.5 is the irregular layout, §6.6 is
  missing data and §6.7 is the sign-flip test. The task brief's "§6.5,
  diagonal-deleted 40 x 40" refers to what is now §6.6.
* **Code under test:** package 0.4.2, the working tree on branch
  `test-suite-paper-fidelity`.

**Conventions adopted (decisions made without the maintainer).**

1. **testthat 3e was added *alongside* the existing base-R suite.** It was not
   converted: rewriting the 16 base-R files would delete files this session
   did not create. The brief makes deletions sign-off-only, and the base-R
   suite is the maintainer's documented contract (`tests/README.md`).
   `testthat (>= 3.0.0)` is in `Suggests`. The reason: it is a test-only
   dependency that the new suite needs. `R CMD check` runs both suites, the
   base-R scripts and `tests/testthat.R`.
2. **Existing assertions are cited rather than duplicated.** The traceability
   table below names both suites.
3. **Tools.** devtools is not installed. The equivalents used are
   `testthat::test_local()` for `devtools::test()`,
   `R CMD check --as-cran` for `devtools::check()`, and
   `roxygen2::roxygenise()` 7.3.2 for `devtools::document()`. covr 3.6.4 and
   roxygen2 7.3.2 come from a local tool library, so nothing was downloaded.
   withr is used only through testthat's own dependency and is never called
   directly.
4. **Branch.** The session works on `test-suite-paper-fidelity`. Its first
   commit, `1b62733`, is a snapshot of the uncommitted 0.4.2 working tree
   exactly as found. `release/0.4.1` is untouched.

---

## Checklist

### Phase 1: Orientation
- [x] Read every file in `R/` (17 files) and `tests/` (16 top level, 8 lower level, 2 helpers, golden)
- [x] Read DESCRIPTION and NAMESPACE
- [x] Read the IPT paper: Procedures 1-2, Algorithms 1-2, Theorems 1 and 4, Proposition 2, Remark 1, Sections 5-6
- [x] Read the RPT paper: Algorithms 1-2 (page images)
- [x] Read the JSS draft: Sections 2-7
- [x] Baseline base-R tests (see below)
- [x] Baseline coverage (see below)
- [x] Baseline `R CMD check --as-cran` (see below)

### Phase 2: Traceability
- [x] Table: paper/draft object -> function:lines -> tests
- [x] Gaps identified

### Phase 3: Tests (testthat 3e, `tests/testthat/`)
- [x] Infrastructure: `tests/testthat.R`, DESCRIPTION `Suggests`/`Config/testthat/edition`
- [x] `helper-reference.R`: naive Procedure 1, seam reconstruction of each front end's group, DGPs, ECDF checker
- [x] `test-procedure1-reference.R`
- [x] `test-algorithm1.R`
- [x] `test-invariances.R`
- [x] `test-pvalue-properties.R`
- [x] `test-confidence-sets.R`
- [x] `test-designs.R`
- [x] `test-bicliques.R`
- [x] `test-dispatch.R` (+ snapshots)
- [x] `test-s3-methods.R` (+ snapshots)
- [x] `test-validation.R`
- [x] `test-reproducibility.R`
- [x] `test-montecarlo.R` (slow, gated)
- [x] `test-draft-replication.R` (+ gravity-gated Sec. 7)

### Phase 4: Comments
- [ ] Header comment on every `R/` file
- [ ] roxygen block (with `@details` step citation) on every function
- [ ] Inline method comments
- [ ] `ARCHITECTURE.md` at root (+ `.Rbuildignore`)
- [ ] Behavior-neutrality proof: `deparse(parse(keep.source = FALSE))` identical per file; roxygenise leaves NAMESPACE unchanged

### Phase 5: Verify and report
- [ ] `test_local()` passes with slow tests skipped
- [ ] ... and with `MWPERM_SLOW_TESTS=true`
- [ ] `R CMD check --as-cran`: no new ERROR/WARNING/NOTE
- [ ] Coverage after
- [ ] Mutation check (a)-(e), throwaway branch deleted
- [ ] Commits: tests and comments separate, small
- [ ] Final report

---

## Phase 1 baseline (2026-09-24, package 0.4.2)

**Base-R suite:** 16/16 top-level files pass (`Rscript tests/test-*.R` after
`R CMD INSTALL .`). `test-lower-level.R` runs the 8 lower-level files. Wall
times in seconds: dyadic 1, equivariance 3, formula 0, golden 9, irregular 0,
layout 1, lower-level 16, main 9, methods 1, missing 1, panel-missing 1,
panel 1, readme 1, signflip 14, threeway 1, validation 1. Total about 60 s.

**Coverage (covr, `type = "tests"`, base-R suite only):** 94.5% overall.

| File | coverable lines | covered | % |
|---|---:|---:|---:|
| R/core.R | 186 | 177 | 95.2 |
| R/dyadic.R | 36 | 36 | 100.0 |
| R/engine.R | 660 | 605 | 91.7 |
| R/formula.R | 44 | 44 | 100.0 |
| R/irregular.R | 177 | 164 | 92.7 |
| R/layout.R | 99 | 99 | 100.0 |
| R/methods.R | 154 | 154 | 100.0 |
| R/missing.R | 373 | 357 | 95.7 |
| R/panel.R | 36 | 36 | 100.0 |
| R/panel_missing.R | 194 | 190 | 97.9 |
| R/perm_set.R | 38 | 38 | 100.0 |
| R/plot.R | 386 | 377 | 97.7 |
| R/signflip.R | 134 | 122 | 91.0 |
| R/threeway.R | 32 | 32 | 100.0 |
| R/unified.R | 688 | 628 | 91.3 |

The uncovered lines that carry the method are:

* the disconnected-set note on the exact path (`engine.R` 394-403);
* the bisection fallback's rejected-centre restart and its island-guard
  widening (`engine.R` 901-904, 954-959);
* the unbounded-bracket return (`engine.R` 918);
* the blockwise `.row_median()` path (`engine.R` 576-585);
* the joint-region budget guards (`engine.R` 1038-1064);
* `.ipt_pvalue()` (`core.R` 433-434);
* the `.cell_code()` 2^53 guard (`core.R` 72-75);
* the exact-biclique budget fallback in `find_bicliques()` (`missing.R` 603);
* the greedy retry-peel success path (`missing.R` 631-638).

**R CMD check --as-cran (baseline):** 0 ERROR, 0 WARNING, 3 NOTEs: "New submission" (CRAN incoming), "unable to verify current time", and HTML-manual tidy warnings from R 4.3.2's Rd2HTML template. The same three NOTEs as the maintainer's documented gate. Check time 118 s. The tarball contains no private or root-level planning files.

---

## Phase 2: Traceability table

Legend. **T:** a testthat file in `tests/testthat/` (the quoted phrase is the
test's name). **B:** a base-R file in `tests/`, with `ll/` for
`lower-level-tests/` and § for the file's own numbered section. Locations are
`file:first line` of the implementing code in the final tree (refreshed
after Phase 4). A row with no test gives its reason.

### IPT paper (Guo, Toulis & Wang 2026)

| # | Paper object | Implementation | Tests |
|---|---|---|---|
| P1 | Eq. (7), stacking cell (i,j) in row (i-1)n+j | front ends build y, D, X (`.make_X` adds the intercept); cells are keyed by `.cell_code` (mixed radix, first coordinate least significant), **not** by lexicographic row order. This is immaterial because the statistic is invariant to a common reordering of the rows. | T: procedure1-reference "the lexicographic stacking of Eq. (7) is immaterial"; B: test-equivariance §1, test-panel §4, test-threeway §3 |
| P2 | Assumption 1 (double exchangeability) | an assumption, not code; stated in `?mwperm_dyadic` | T: montecarlo t(3), Cauchy, Moulton (size under GTW Eq. 8 errors) |
| P3 | Eq. (9): G = {(pi_k, sigma_k)}, paired by k | `mwperm_dyadic()` perm_builder pairs `Grow[[k]]` with `Gcol[[k]]`; `.build_obs_perms()` | T: procedure1-reference (all), designs; B: ll/test-obsperms §1 |
| P4 | y_{pi,sigma}: entry (i,j) is y_{pi(i) sigma(j)} | `.build_obs_perms()` (gather = `match()` or position table on the mapped cell code) | T: procedure1-reference "the package builds y_{pi_k, sigma_k} as GTW define it" |
| P5 | Proc. 1 step 1: V_k is an orthonormal basis of col([X, X_k])^perp | `.ipt_prepare()`: residualizes D on [X, X_k] with an eigen pseudo-inverse of the 2p x 2p Gram matrix, cut at 1e-14 (the rank, not 2p) | T: procedure1-reference (a_k, b_k to 1e-8; row/col FE rank deficiency), designs (panel time FE); B: ll/test-projection §1-4 |
| P6 | a_k = norm(D' V_k V_k' y), b_k = norm(D' V_k V_k' y_k) | `.ipt_prepare()` caches u = Dr'y, v = Dr'y_k, M = Dr'Dr, W = Dr'D_k; `.ipt_eval()` evaluates a = abs(u - M b), b = abs(v - W b), or the norms when d > 1 | T: procedure1-reference; B: ll/test-projection §3 |
| P7 | Eq. (10): minorized p-value, min over j = 1..K, indicator <= | `.ipt_eval()` `(1 + sum(b >= amin)) / Kp1`; `.pval_matrix()` (vectorised) | T: procedure1-reference (identical p-values; "min over the K non-identity elements"), pvalue-properties (grid, ties); B: ll/test-pvalue §1-5 |
| P8 | testing beta = b means running Proc. 1 on y - D b | the affine form in `.ipt_eval()`; `beta0` in `.ipt_engine()` | T: procedure1-reference "a non-zero null ...", invariances "testing beta = b ..."; B: test-equivariance §3 |
| P9 | Proc. 1 step 3: CI = {b : pval(b) > alpha} | `.invert_ci()` (exact via `.ci_breakpoints()` + `.exact_ci_set()`; grid; bisection), `.invert_region()` when d > 1 | T: confidence-sets (duality at every end point, grid and bisection versus exact, region); **D10**; B: ll/test-exact-ci, ll/test-invert-ci-grid |
| P10 | Theorem 1: P(pval <= alpha given X, D) <= alpha, for p < N/2 | engine guard `N <= 2p` stops | T: validation "p >= N/2 is refused", montecarlo (ECDF at every atom) |
| P11 | Remark 1: median over repetitions | `.agg_pvals()`, `.row_median()`; `pvalue` in `.ipt_engine()` | T: pvalue-properties (median, odd and even; median2); B: ll/test-aggregate |
| P12 | Algorithm 1 | `build_perm_set()` | T: algorithm1 (all); **D1, D2**; B: ll/test-permset |
| P13 | Proposition 2 (closure) | `build_perm_set()`: the elements are the powers of one generator | T: algorithm1 "closed under composition ...", "psi_k is psi_1 composed ..."; B: ll/test-permset §2, ll/test-obsperms (observation level) |
| P14 | Definition 1 (fully observed, disjoint blocks) | `find_bicliques()` | T: bicliques "blocks are fully observed, disjoint ..."; B: test-missing §1 |
| P15 | Assumption 4 (mask independent of errors) | an assumption; `.choose_common_levels()` reads the mask only | T: montecarlo MCAR |
| P16 | Procedure 2, steps 1-4 | `mwperm_missing()`, `.build_obs_perms_blocks()` | T: designs "missing: Procedure 2 on the pooled blocks ...", "... never mixed"; B: test-missing, ll/test-obsperms §5 |
| P17 | Theorem 4 | -- | T: montecarlo MCAR (ECDF, conditional on one mask) |
| P18 | Algorithm 2 (App. A) | `find_bicliques()` (`.grow_biclique()` greedy default, `.max_biclique_exact()` branch and bound); documented deviations below | T: bicliques (brute-force maximum, validity, fallback) |
| P19 | Sec. 6.1, InvA, three-way | `mwperm_threeway()` | T: designs three-way; B: test-threeway |
| P20 | Sec. 6.2, InvB, panel | `mwperm_panel()` (time's group is NULL, i.e. held fixed) | T: designs panel (x2), invariances trend, montecarlo panel and negative control; B: test-panel |
| P21 | Sec. 6.3, layout | `mwperm_layout()`, `.build_obs_perms_layout()`, `.within_cell_slot()`, `.downsample_to_L0()` | T: designs layout (x3); **D7**; B: test-layout |
| P22 | Sec. 6.4, irregular | `mwperm_irregular()`, `.irregular_design()` | T: designs irregular, draft-replication Sec. 6.5; B: test-irregular |
| P23 | Sec. 9 (open): incomplete panels | `mwperm_panel_missing()`, `.panel_missing_design()`, `.choose_common_levels()` | T: designs incomplete panel; B: test-panel-missing |
| P24 | revised Assumption 2 (double sign symmetry) and its Section E procedure | `mwperm_dyadic_het()`, `build_flip_set()`, `.build_obs_flips()`, `.apply_op()` | B: test-signflip (against the authors' port); T: pvalue-properties (grid 2^(n_flip-1)), dispatch, reproducibility. **Gap:** Section E is not available, so fidelity to the paper's own procedure is unchecked |
| P25 | Theorems 2, 3, 5 and Proposition 1 (power; asymptotic) | -- | **No test, by nature:** asymptotic statements with no finite-sample assertion. T: montecarlo "power ..." is a sanity check only |

### RPT paper (Wen, Wang & Wang 2025)

| # | Paper object | Implementation | Tests |
|---|---|---|---|
| R1 | Algorithm 1, the permutation set, with a repeat-until loop on tr(V_0 V_0' P_k) | `build_perm_set()` follows GTW's Algorithm 1, which drops the loop (one draw of pi) | T: algorithm1; documented deviation (GTW's choice, not the package's) |
| R2 | Algorithm 2, RPT: V_k = P_k V_0, V~_k in span(V_0) intersect span(V_k) | `.ipt_prepare()`: the same subspace, col([X, X_k])^perp | T: procedure1-reference |

### JSS draft (2026-09-18 build)

| # | Draft statement | Implementation | Tests |
|---|---|---|---|
| J1 | Sec. 2.2, Eq. vk: rank form; eigen cut 1e-14 | `.ipt_prepare()` | T: procedure1-reference rank-deficient case, designs panel time FE |
| J2 | Sec. 2.3, Eq. aggset: one aggregation rule for the p-value and the set | `.agg_pvals()`, used in `.ipt_engine()` and on every inversion path | T: confidence-sets "the reported p-value is the one the set inverts"; B: ll/test-aggregate |
| J3 | Sec. 2.4, Eq. breaks; closure convention | `.ci_breakpoints()`, `.exact_ci_set()` | T: confidence-sets (end points are roots; duality); **D10** |
| J4 | Sec. 2.5: resolution 1/(K+1) <= alpha; median2 doubles the floor; default K = min - 1, capped at 199; an even-rep median can leave the grid | `.default_K()`, `p_floor` in `.ipt_engine()`, the notes | T: confidence-sets "resolution guard", pvalue-properties "median ..."; B: test-validation §5 |
| J5 | Sec. 2.7: extensions | front ends | T: designs |
| J6 | Sec. 3.2, Table 1: dispatch | `mwperm_check()`, `mwperm()` | T: dispatch (one test per row); B: test-main §4 |
| J7 | Sec. 3.2, Table 2: negative control | -- | T: montecarlo "negative control" |
| J8 | Sec. 3.2: forcing threeway warns on a time-like index | `mwperm_check()` | B: test-main §3c |
| J9 | Sec. 3.3: `mwperm_check()` printed output | `mwperm_check()`, `print.mwperm_design()` | T: draft-replication Sec. 3.3, dispatch snapshots |
| J10 | Sec. 3.4: first example; `mwperm()` identical to the worker | `mwperm()`, `mwperm_dyadic()`, `print.mwperm()` | T: draft-replication Sec. 3.4 |
| J11 | Sec. 3.5: slots and the six methods | the `structure()` in `.ipt_engine()`, methods.R, formula.R | T: s3-methods; B: test-methods |
| J12 | Sec. 3.6: `build_perm_set()`/`build_flip_set()` printed output; the prose | `build_perm_set()`, `build_flip_set()` | T: draft-replication Sec. 3.6, algorithm1; **D3** |
| J13 | Sec. 3.7: conf_set, ci_method, empty set, root budget, `confint()` error | `.ipt_engine()`, `.invert_ci()`, `confint.mwperm()` | T: confidence-sets, draft-replication Sec. 3.7; **D6**; B: test-validation §5-6 |
| J14 | Sec. 3.8: seed + r - 1; rep_seed x stride + offset; RNG restored; parallel identical to serial | `.ipt_engine()` (`seeds`), `.sub_seed()`, `.plapply()`, `.save_seed()`/`.restore_seed()` | T: reproducibility |
| J15 | Sec. 3.9: biclique finder | `find_bicliques()` | T: bicliques |
| J16 | Sec. 4: the synthetic data sets | data/, data-raw/ | T: draft-replication Sec. 4 |
| J17 | Sec. 5.1: Moulton design (Table 3) | -- | T: montecarlo Moulton. **Table 3's numbers are not asserted:** they are not reproducible from the printed snippet (prior audit, finding P3) |
| J18 | Sec. 5.2: verification record | -- | T: procedure1-reference, designs, bicliques cover the same ground |
| J19 | Secs. 6.1-6.7: printed outputs | front ends | T: draft-replication (all pass numerically); **D4, D5, D6** |
| J20 | Sec. 7: real data | `mwperm()` on the gravity data | T: draft-replication Sec. 7 (runs only if gravity is installed); **Sec. 7.4 (CEPII) is skipped:** the data are not distributed |

### Documented deviations (the code differs from the printed text on purpose; not discrepancies)

* **Rank instead of N - 2p** (GTW Procedure 1 step 1; draft Sec. 2.2). The
  stack [X | X_k] is rank-deficient by construction: the permutation fixes the
  intercept, and panel time dummies too. The package works with its numerical
  rank. The reference test agrees to 1e-8.
* **Algorithm 2's stopping rule.** GTW loop "while M has nonzero entries" and
  solve for the maximum biclique. `find_bicliques()` stops at the first block
  below `min_block`, and its default is greedy. Remark 2 of GTW makes an
  approximate maximum acceptable, so this costs power only.
* **WWW's repeat-until loop in Algorithm 1** is dropped, as GTW's Algorithm 1
  drops it.
* **Package extensions beyond the paper:**
  * `time_fe = TRUE` (period dummies are invariant under InvB);
  * `aggregate = "median2"`;
  * `permute = "rows"/"cols"` (a subgroup);
  * `L0` in `mwperm_layout()` and `mwperm_panel_missing()`;
  * the incomplete-panel design;
  * incomplete arrays for the sign-flip test.

  Each is documented in its help page.

### Gaps (with reasons)

* P25: asymptotic power theorems. There is no finite-sample assertion to test.
* P24: the revised paper's Section E is not available, so `build_flip_set()`
  cannot be checked against the paper's own construction. The base suite
  checks it against the authors' script instead.
* The disconnected-set note on the exact path, and the bisection fallback's
  island guard, have no front-end fixture: no seeded fit out of about 1,200
  tried produced a disconnected set. Base-suite internals cover the helpers,
  not the note text. Seam F-3 below would close this.
* J17: the Table 3 numbers (see above).
* J20, Sec. 7.4: the CEPII data are unavailable.

---

## Discrepancy log

Each entry gives the inputs, the expected and observed values, the reference,
and a verdict. Skipped tests carry the ID, and `MWPERM_SHOW_DISCREPANCIES=true`
runs them.

**D1: Algorithm 1 prints "i mod (K+1)".**
* Reference: GTW Algorithm 1 display; WWW Algorithm 1, Eq. (10), has the
  same wording.
* Expected (the formula as printed, with the 0-based R residue): psi~_k is a
  permutation of 1..n.
* Observed (n = 6, K = 2, k = 1): the images are `2 3 4 5 6 7`. Index K+1
  of each block is sent into the next block, and index n is sent to n+1.
  The code uses the 1-based residue ((i-1) mod (K+1)) + 1 (`pos0` in
  `build_perm_set()`), and element k equals pi^{-1} o psi~_k o pi exactly
  under that reading (T: algorithm1).
* Verdict: a typo in the paper. **The paper should change**
  (write ((i-1) mod (K+1)) + 1). The code is correct.
* Test: algorithm1 "Algorithm 1 as printed ..." (skipped).

**D2: Algorithm 1's worked example lists the inverse map.**
* Reference: GTW p. 11-12; WWW p. 736.
* The example gives S_1^k = (K+2-k, ..., K+1, 1, ..., K+1-k) as
  (psi~_k(1), ..., psi~_k(K+1)). The displayed formula gives psi~_k(1) = 1+k.
  At n = 8, K = 3, k = 1 the formula gives `2 3 4 1` and the example `4 1 2 3`.
* The two differ by k <-> K+1-k, i.e. by inversion, so the SET, which is
  what validity uses, is identical (tested). The code follows the formula.
* Verdict: **the paper should align the example with the formula.** No code
  change.
* Test: algorithm1 "worked example ..." (skipped).

**D3: Draft Sec. 3.6 prose about `build_perm_set(n = 6, K = 2, seed = 1)`.**
* The draft says the elements shift "the two blocks {1,2,3} and {4,5,6} (here
  under the identity relabeling drawn from seed = 1)".
* Observed:
  * for seed 1 the relabelling is pi = (1, 4, 3, 6, 2, 5), not the identity;
  * the printed element `5 6 1 2 3 4`, which the code reproduces exactly,
    has orbits {1,3,5} and {2,4,6}, which are pi^{-1}({1,2,3}) and
    pi^{-1}({4,5,6}).
* `attr(, "block_size")` is K+1 = 3. That is the block length in
  *relabelled* coordinates, and also the group order.
* Verdict: **the draft should change the sentence.** Suggested wording:
  "cyclically shifts the blocks {1,3,5} and {2,4,6}, the pre-images of
  {1,2,3} and {4,5,6} under the relabelling drawn from seed = 1". The code
  and the printed output are right.
* Test: algorithm1 "draft Sec. 3.6 prose ..." (skipped).

**D4: Draft Sec. 6.7, printed header of a sign-flip fit.**
* The draft prints "Invariant permutation test (mwperm)". Since 0.4.2 the
  package prints "Invariant sign-flip test (mwperm)", a deliberate change
  recorded in NEWS 0.4.2. Every number in the block reproduces.
* Verdict: **the draft should change.**
* Test: draft-replication "Sec. 6.7: ... printed header" (skipped).

**D5: Draft Secs. 2.7, 3.2, 3.3 (Table 2 caption) and 6.5: `trim = "levels"` and the 0.4.1 irregular note.**
* `mwperm_irregular(trim = )` was removed in 0.4.2, and the first note was
  rewritten to say the random trim is exact only for exchangeable
  replicates. This follows the advisors' decision of 2026-09-22 and the
  2026-09-23 audit.
* Every number in the Sec. 6.5 block reproduces.
* Verdict: **the draft should change.** Describe
  `mwperm_panel_missing(time =, L0 =)` for periods, and reprint the note.
* Tests: draft-replication "Sec. 6.5: ... printed note" and
  "... trim = \"levels\" option" (both skipped).

**D6: Draft Secs. 2.7, 3.7 and 6.7 say the sign-flip default is n_flip = 8.**
* 0.4.2 made the default the resolution rule: 6 at alpha = 0.05 (7 under
  median2).
* The draft's Sec. 3.7 budget claim (n_flip = 8, 10 reps: bisection) still
  holds for an explicit `n_flip = 8` (tested).
* Verdict: **the draft should change.**
* Test: draft-replication "the default n_flip is 8" (skipped).

**D7: GTW Sec. 6.3 claims layout validity under a shared i.i.d. zeta_l.**
* GTW, and the draft's Sec. 2.7, say the layout test is valid for
  eps_ijl = eta_ij + zeta_l + u_ijl with zeta_l i.i.d. and shared across
  cells. The draft's Sec. 6.4 example DGP includes exactly such a
  `zeta[idl]`.
* Step (i) applies Algorithm 1 to each cell independently. So group element
  k permutes l differently in different cells, and the array
  (zeta_{pi^c_k(l)}) is not equal in law to (zeta_l): the invariance
  Theorem 1 needs fails. Observed (3 x 3 cells of 6, K = 5): the five
  non-identity elements apply 9 different within-cell maps.
* `?mwperm_layout` already says this component "is not covered" and that
  simulations found no measurable size effect.
* Verdict: **the paper and the draft should change.** Either drop the
  zeta_l example, or state that it needs one common within-cell permutation
  (possible when cell sizes are equal). Validity under the paper's other
  example (eta_ij + u_ijl) is unaffected. No code change; a Monte Carlo that
  quantifies the size effect is follow-up F-6.
* Test: designs "layout under a shared replicate effect zeta_l" (skipped).

**D10: Reported CI end points are closures, not members of {b : pval(b) > alpha}.**
* Reference: GTW Procedure 1 step 3 defines CI = {b : pval(b) > alpha}. The
  exact path reports each connected component closed, so an end point can be
  a rejected jump.
* Observed: for the `trade_dyadic` anchor (`mwperm_dyadic`, seed 1,
  n_reps 10), the lower end is -1.251351 with p = 0.05 = alpha, i.e.
  rejected, while p(-1.251351 + 1e-12) = 0.0625.
* The direction is outward: the reported set is never smaller than the true
  one, and just inside and just outside every end point behave correctly
  (tested). The package documents the convention (`?confint.mwperm`, NEWS
  0.3.0), and so does the draft (Sec. 2.4).
* Verdict: **documented and conservative, so keep the code.** A change would
  report the attained side, which moves seeded end points: a numeric change
  that needs sign-off.
* Test: confidence-sets "reported end points are in ..." (skipped).

(D8 and D9 are unused. The Algorithm 2 stopping rule and the even-rep median
turned out to be documented behavior, not discrepancies. See *Documented
deviations*, and the draft's Sec. 2.5, which already says the even-rep
median can leave the grid.)

---

## Follow-ups (testing seams, not implemented)

* **F-1.** Expose each front end's per-repetition group. Examples: an
  internal `.dyadic_design()` like the existing `.panel_missing_design()`
  and `.irregular_design()`, or a `keep_groups = TRUE` debug field. The
  tests currently rebuild the groups from the documented seed scheme
  (`seed + r - 1`, `rep_seed * stride + j`), which duplicates that scheme in
  `helper-reference.R`.
* **F-2.** Optionally keep the per-repetition cross products (`u, v, M, W`)
  on the fit, e.g. `keep_prep = TRUE`. The reference comparison could then
  use the fit itself, not `.ipt_prepare()` on a rebuilt group. This would
  also enable the reserved `plot(type = "null"/"profile")`.
* **F-3.** An internal entry that runs the engine's inversion and note logic
  on a supplied `prep_list`, so that the disconnected-set note and the
  bisection island guard can be exercised end to end.
* **F-4.** Record the rows kept by `mwperm_layout(L0 = )` on the fit, so the
  balancing draw does not have to be re-derived through
  `.downsample_to_L0()`.
* **F-5.** A seeded Monte Carlo of the layout test under a shared zeta_l
  (D7), to quantify what the documentation now says qualitatively.
* **F-6.** Replicate draft Sec. 7.4 once the CEPII database can be supplied
  locally.
