# TESTING_PLAN.md: paper-to-code traceability and discrepancy log

This file records which code implements each object of Guo, Toulis & Wang
(2026) and Wen, Wang & Wang (2025), and which tests check it. It also logs
every place where the code differs from those papers or from the JSS draft,
whether on purpose or not. It is listed in `.Rbuildignore`, so it is not in
the tarball. The name is historical: the file began as the working plan of
the session that built the testthat suite, and that session's record is kept
as an appendix.

| Section | Use it to |
|---|---|
| [Traceability table](#traceability-table) | find the function and the tests behind a paper object (P: IPT paper, R: RPT paper, J: JSS draft) |
| [Documented deviations](#documented-deviations-the-code-differs-from-the-printed-text-on-purpose-not-discrepancies) | see where the code departs from the printed text on purpose, and why |
| [Gaps](#gaps-with-reasons) | see which paper objects have no test, and why |
| [Discrepancy log](#discrepancy-log) | look up a `Discrepancy D<n>` skip: inputs, expected and observed values, verdict |
| [Follow-ups](#follow-ups-testing-seams-not-implemented) | testing seams worth adding (F-1 to F-6) |
| [Appendix](#appendix-the-session-that-built-the-suite-2026-09-24) | the 2026-09-24 session: its conventions, baseline and final report, including the decisions awaiting the maintainer's review |

How the two test suites are organised, and how to run them, is in
`tests/README.md`.

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

---

## Traceability table

Legend. **T:** a testthat file in `tests/testthat/` (the quoted phrase is the
test's name). **B:** a base-R file in `tests/`, with `ll/` for
`lower-level-tests/` and § for the file's own numbered section. Locations are
`file:line` in the final tree, commit `5ef3343` (after the comment pass). A row with no test gives its reason.

### IPT paper (Guo, Toulis & Wang 2026)

| # | Paper object | Implementation | Tests |
|---|---|---|---|
| P1 | Eq. (7), stacking cell (i,j) in row (i-1)n+j | front ends build y, D, X (`.make_X`, engine.R:1230, adds the intercept); cells are keyed by `.cell_code` (core.R:101) (mixed radix, first coordinate least significant), **not** by lexicographic row order. This is immaterial because the statistic is invariant to a common reordering of the rows. | T: procedure1-reference "the lexicographic stacking of Eq. (7) is immaterial"; B: test-equivariance §1, test-panel §4, test-threeway §3 |
| P2 | Assumption 1 (double exchangeability) | an assumption, not code; stated in `?mwperm_dyadic` | T: montecarlo t(3), Cauchy, Moulton (size under GTW Eq. 8 errors) |
| P3 | Eq. (9): G = {(pi_k, sigma_k)}, paired by k | `mwperm_dyadic()` perm_builder (dyadic.R:255) pairs `Grow[[k]]` with `Gcol[[k]]`; `.build_obs_perms()` (core.R:570) | T: procedure1-reference (all), designs; B: ll/test-obsperms §1 |
| P4 | y_{pi,sigma}: entry (i,j) is y_{pi(i) sigma(j)} | `.build_obs_perms()` (core.R:570; gather = `match()` or position table on the mapped cell code) | T: procedure1-reference "the package builds y_{pi_k, sigma_k} as GTW define it" |
| P5 | Proc. 1 step 1: V_k is an orthonormal basis of col([X, X_k])^perp | `.ipt_prepare()` (core.R:297; the eigen pseudo-inverse at core.R:377): residualizes D on [X, X_k], cut at 1e-14 (the rank, not 2p) | T: procedure1-reference (a_k, b_k to 1e-8; row/col FE rank deficiency), designs (panel time FE); B: ll/test-projection §1-4 |
| P6 | a_k = norm(D' V_k V_k' y), b_k = norm(D' V_k V_k' y_k) | `.ipt_prepare()` (core.R:297) caches u = Dr'y, v = Dr'y_k, M = Dr'Dr, W = Dr'D_k; `.ipt_eval()` (core.R:465) evaluates a = abs(u - M b), b = abs(v - W b), or the norms when d > 1 | T: procedure1-reference; B: ll/test-projection §3 |
| P7 | Eq. (10): minorized p-value, min over j = 1..K, indicator <= | `.ipt_eval()` `(1 + sum(b >= amin)) / Kp1` (core.R:494); `.pval_matrix()` (engine.R:661, vectorised) | T: procedure1-reference (identical p-values; "min over the K non-identity elements"), pvalue-properties (grid, ties); B: ll/test-pvalue §1-5 |
| P8 | testing beta = b means running Proc. 1 on y - D b | the affine form in `.ipt_eval()` (core.R:465); `beta0` in `.ipt_engine()` (engine.R:204) | T: procedure1-reference "a non-zero null ...", invariances "testing beta = b ..."; B: test-equivariance §3 |
| P9 | Proc. 1 step 3: CI = {b : pval(b) > alpha} | `.invert_ci()` (engine.R:873; exact via `.ci_breakpoints()` engine.R:719 + `.exact_ci_set()` engine.R:775; grid; bisection), `.invert_region()` (engine.R:1063) when d > 1 | T: confidence-sets (duality at every end point, grid and bisection versus exact, region); **D10**; B: ll/test-exact-ci, ll/test-invert-ci-grid |
| P10 | Theorem 1: P(pval <= alpha given X, D) <= alpha, for p < N/2 | engine guard `N <= 2p` stops (engine.R:261) | T: data-validation "p >= N/2 is refused", montecarlo (ECDF at every atom) |
| P11 | Remark 1: median over repetitions | `.agg_pvals()` (engine.R:561), `.row_median()`; `pvalue` (engine.R:381) | T: pvalue-properties (median, odd and even; median2); B: ll/test-aggregate |
| P12 | Algorithm 1 | `build_perm_set()` (perm_set.R:93; psi_k at perm_set.R:148) | T: algorithm1 (all); **D1, D2**; B: ll/test-permset |
| P13 | Proposition 2 (closure) | `build_perm_set()`: the elements are the powers of one generator | T: algorithm1 "closed under composition ...", "psi_k is psi_1 composed ..."; B: ll/test-permset §2, ll/test-obsperms (observation level) |
| P14 | Definition 1 (fully observed, disjoint blocks) | `find_bicliques()` (missing.R:594) | T: bicliques "blocks are fully observed, disjoint ..."; B: test-missing §1 |
| P15 | Assumption 4 (mask independent of errors) | an assumption; `.choose_common_levels()` reads the mask only | T: montecarlo MCAR |
| P16 | Procedure 2, steps 1-4 | `mwperm_missing()` (missing.R:152), `.build_obs_perms_blocks()` (missing.R:352) | T: designs "missing: Procedure 2 on the pooled blocks ...", "... never mixed"; B: test-missing, ll/test-obsperms §5 |
| P17 | Theorem 4 | -- | T: montecarlo MCAR (ECDF, conditional on one mask) |
| P18 | Algorithm 2 (App. A) | `find_bicliques()` (missing.R:594; `.grow_biclique()` missing.R:723 greedy default, `.max_biclique_exact()` missing.R:849 branch and bound); documented deviations below | T: bicliques (brute-force maximum, validity, fallback) |
| P19 | Sec. 6.1, InvA, three-way | `mwperm_threeway()` (threeway.R:91) | T: designs three-way; B: test-threeway |
| P20 | Sec. 6.2, InvB, panel | `mwperm_panel()` (panel.R:120; perm_builder panel.R:164, time's group is NULL, i.e. held fixed) | T: designs panel (x2), invariances trend, montecarlo panel and negative control; B: test-panel |
| P21 | Sec. 6.3, layout | `mwperm_layout()` (layout.R:146), `.build_obs_perms_layout()` (layout.R:290), `.within_cell_slot()` (core.R:167), `.downsample_to_L0()` (layout.R:265) | T: designs layout (x3); **D7**; B: test-layout |
| P22 | Sec. 6.4, irregular | `mwperm_irregular()` (irregular.R:192), `.irregular_design()` (irregular.R:363) | T: designs irregular, draft-replication Sec. 6.5; B: test-irregular |
| P23 | Sec. 9 (open): incomplete panels | `mwperm_panel_missing()` (panel_missing.R:186), `.panel_missing_design()` (panel_missing.R:387), `.choose_common_levels()` (panel_missing.R:336) | T: designs incomplete panel; B: test-panel-missing |
| P24 | revised Assumption 2 (double sign symmetry) and its Section E procedure | `mwperm_dyadic_het()` (signflip.R:523), `build_flip_set()` (signflip.R:170), `.build_obs_flips()` (signflip.R:318), `.apply_op()` (core.R:217) | B: test-signflip (against the authors' port); T: pvalue-properties (grid 2^(n_flip-1)), dispatch, reproducibility. **Gap:** Section E is not available, so fidelity to the paper's own procedure is unchecked |
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
| J4 | Sec. 2.5: resolution 1/(K+1) <= alpha; median2 doubles the floor; default K = min - 1, capped at 199; an even-rep median can leave the grid | `.default_K()` (engine.R:1262), `p_floor` in `.ipt_engine()`, the notes | T: confidence-sets "resolution guard", pvalue-properties "median ..."; B: test-validation §5 |
| J5 | Sec. 2.7: extensions | front ends | T: designs |
| J6 | Sec. 3.2, Table 1: dispatch | `mwperm_check()` (unified.R:240), `mwperm()` (unified.R:1046) | T: dispatch (one test per row); B: test-main §4 |
| J7 | Sec. 3.2, Table 2: negative control | -- | T: montecarlo "negative control" |
| J8 | Sec. 3.2: forcing threeway warns on a time-like index | `mwperm_check()` | B: test-main §3c |
| J9 | Sec. 3.3: `mwperm_check()` printed output | `mwperm_check()`, `print.mwperm_design()` | T: draft-replication Sec. 3.3, dispatch snapshots |
| J10 | Sec. 3.4: first example; `mwperm()` identical to the worker | `mwperm()`, `mwperm_dyadic()`, `print.mwperm()` | T: draft-replication Sec. 3.4 |
| J11 | Sec. 3.5: slots and the six methods | the `structure()` in `.ipt_engine()`, methods.R, formula.R | T: s3-methods; B: test-methods |
| J12 | Sec. 3.6: `build_perm_set()`/`build_flip_set()` printed output; the prose | `build_perm_set()`, `build_flip_set()` | T: draft-replication Sec. 3.6, algorithm1; **D3** |
| J13 | Sec. 3.7: conf_set, ci_method, empty set, root budget, `confint()` error | `.ipt_engine()`, `.invert_ci()`, `confint.mwperm()` | T: confidence-sets, draft-replication Sec. 3.7; **D6**; B: test-validation §5-6 |
| J14 | Sec. 3.8: seed + r - 1; rep_seed x stride + offset; RNG restored; parallel identical to serial | `.ipt_engine()` (`seeds`, engine.R:330), `.sub_seed()` (engine.R:1436), `.plapply()` (engine.R:56), `.save_seed()`/`.restore_seed()` (perm_set.R:165) | T: reproducibility |
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
  `helper-reference.R`. The mutation check shows the cost. A
  fresh (pi, sigma) per panel period (mutation d) was caught only by the two
  p-value comparisons, because the structural test "one (pi, sigma) in every
  period" had to build its own group.
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

---

## Appendix: the session that built the suite (2026-09-24)

A record of how the testthat suite and the source-comment pass were
produced. The numbers here describe the tree as it was that day; the
sections above are the part kept current.

### Conventions adopted (decisions made without the maintainer)

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

### Checklist

#### Phase 1: Orientation
- [x] Read every file in `R/` (17 files) and `tests/` (16 top level, 8 lower level, 2 helpers, golden)
- [x] Read DESCRIPTION and NAMESPACE
- [x] Read the IPT paper: Procedures 1-2, Algorithms 1-2, Theorems 1 and 4, Proposition 2, Remark 1, Sections 5-6
- [x] Read the RPT paper: Algorithms 1-2 (page images)
- [x] Read the JSS draft: Sections 2-7
- [x] Baseline base-R tests (see below)
- [x] Baseline coverage (see below)
- [x] Baseline `R CMD check --as-cran` (see below)

#### Phase 2: Traceability
- [x] Table: paper/draft object -> function:lines -> tests
- [x] Gaps identified

#### Phase 3: Tests (testthat 3e, `tests/testthat/`)
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
- [x] `test-data-validation.R` (named `test-validation.R` until 2026-09-24; renamed so it no longer shares a name with the base-R file)
- [x] `test-reproducibility.R`
- [x] `test-montecarlo.R` (slow, gated)
- [x] `test-draft-replication.R` (+ gravity-gated Sec. 7)

#### Phase 4: Comments
- [x] Header comment on every `R/` file
- [x] roxygen block (with `@details` step citation) on every function
- [x] Inline method comments
- [x] `ARCHITECTURE.md` at root (+ `.Rbuildignore`)
- [x] Behavior-neutrality proof: `deparse(parse(keep.source = FALSE))` identical per file; roxygenise leaves NAMESPACE unchanged

#### Phase 5: Verify and report
- [x] `test_local()` passes with slow tests skipped
- [x] ... and with `MWPERM_SLOW_TESTS=true`
- [x] `R CMD check --as-cran`: no new ERROR/WARNING/NOTE
- [x] Coverage after
- [x] Mutation check (a)-(e), throwaway branch deleted
- [x] Commits: tests and comments separate, small
- [x] Final report

### Phase 1 baseline (2026-09-24, package 0.4.2)

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

### Final report (2026-09-24)

Written for a reader who saw none of the work. The evidence (scripts, logs,
raw tables) is cached locally under the session's simulation-cache
directory, so every number can be re-derived without re-running anything.

#### What was done

* **Test suite.** A testthat 3e suite, `tests/testthat/`: 13 test files, one
  helper, two snapshot files, 111 tests. It was added *alongside* the
  existing base-R suite (16 top-level scripts plus 8 lower-level ones),
  which is untouched and still runs.
  * The heart of it is `helper-reference.R`: a deliberately naive
    Procedure 1 (a QR residual maker with numerical rank, then a_k, b_k
    and Eq. 10), written from the paper.
  * The package's statistics are compared with it on the same group, and
    so is every per-repetition p-value, the group being rebuilt from the
    seed. This covers the dyadic, three-way, panel, layout, missing,
    incomplete-panel and irregular designs, d > 1, non-zero nulls, and
    rank-deficient [X | X_k] (intercept, row and column fixed effects,
    time fixed effects).
* **Comments.**
  * Every `R/` file opens with a header: its purpose, its paper sections,
    and its pipeline stage.
  * All 84 functions carry an `@details` line citing the paper step they
    implement, or saying they are not one.
  * Inline comments name the paper symbol on every line that carries the
    method.
  * `ARCHITECTURE.md` gives the pipeline, a notation table and a
    file-to-paper map.
  * The pass fixed four stale comments and one wrong help sentence
    (`?mwperm_check`, on `design = "dyadic_het"`).
* **Behavior-neutrality.** Both checks hold:
  * `deparse(parse(f, keep.source = FALSE))` is `identical()` before and
    after for all 17 files;
  * roxygen2 7.3.2 regenerates a `NAMESPACE` identical to the shipped one,
    both before and after.

  The only `man/` change is the new `@details` paragraph on the exported
  pages (hand-applied exactly as roxygen appends it; `checkRd()` clean),
  plus the one corrected sentence.

#### Results at a glance

| Check | Result |
|---|---|
| base-R suite (16 top-level files) | 16/16 pass, before and after (about 60 s) |
| testthat, default (`NOT_CRAN=true`) | 111 tests, 0 failed, 17 skipped (7 slow, 9 discrepancy, 1 CEPII), about 30 s |
| testthat, `MWPERM_SLOW_TESTS=true` | 111 tests, 0 failed, 10 skipped (9 discrepancy, 1 CEPII), 185 s wall with other jobs running (the Monte Carlo alone: 140 s) |
| `R CMD check --as-cran` | 0 ERROR, 0 WARNING, 3 NOTEs, the same three as the baseline (new submission; clock; HTML-manual tidy). `testthat.R` runs in 26 s under check |
| coverage (covr, both suites under check conditions) | 94.5% -> 95.3% |

Monte Carlo (Theorem 1). One repetition per test, which is what the theorem
covers; the empirical CDF is checked at every atom j/(K+1) against
atom + 3 SE. Every design passes at every atom.

| Design | Sims | Size at 0.05 (or at the 0.04 atom) |
|---|---:|---|
| dyadic, two-way random effects, t_3 errors | 1000 | 0.033 |
| dyadic, two-way random effects, Cauchy errors | 1000 | 0.035 |
| draft Sec. 5.1 Moulton design (node-level d) | 1000 | 0.031 (atom 0.04; naive OLS rejects > 20% there, as asserted) |
| MCAR mask, Procedure 2, K = 19 | 1000 | 0.044 |
| panel with a random-walk trend and AR(1) errors | 1000 | 0.045 |

Negative control (draft Sec. 3.2), 200 simulations at alpha = 0.05: the
three-way test rejects 0.690 and the panel test 0.060 (its bound is
0.096). Power at beta = 0.6 is 1.000.

#### Coverage per file (covr lines; before = base-R suite only, after = both suites)

| File | before | after |
|---|---:|---:|
| R/core.R | 95.2 | 95.2 |
| R/dyadic.R | 100.0 | 100.0 |
| R/engine.R | 91.7 | 93.0 |
| R/formula.R | 100.0 | 100.0 |
| R/irregular.R | 92.7 | 92.7 |
| R/layout.R | 100.0 | 100.0 |
| R/methods.R | 100.0 | 100.0 |
| R/missing.R | 95.7 | 97.3 |
| R/panel.R | 100.0 | 100.0 |
| R/panel_missing.R | 97.9 | 97.9 |
| R/perm_set.R | 100.0 | 100.0 |
| R/plot.R | 97.7 | 97.7 |
| R/signflip.R | 91.0 | 91.0 |
| R/threeway.R | 100.0 | 100.0 |
| R/unified.R | 91.3 | 92.9 |
| **total** | **94.5** | **95.3** |

The "after" figure runs under R CMD check conditions, where the slow and
snapshot tests are skipped. The coverable-line counts are identical before
and after, as they must be if the comment pass changed no code.

#### Traceability gaps

See *Gaps (with reasons)* above. In short:
* asymptotic power theorems (untestable in finite samples);
* the revised paper's Section E (not available);
* the disconnected-set note and the island guard (no natural fixture;
  seam F-3);
* the draft's Table 3 (not reproducible from its snippet);
* the CEPII panel (data not distributed).

#### Discrepancies

None is a code bug, and no finding changes a number. Details and evidence
are in the *Discrepancy log* above.

| ID | Where | Who should change |
|---|---|---|
| D1 | GTW Alg. 1 display, "i mod (K+1)" (0-based reading is not a permutation) | paper |
| D2 | GTW Alg. 1 worked example lists the inverse map (same set) | paper |
| D3 | draft Sec. 3.6 prose (the seed-1 relabelling is not the identity; the orbits are {1,3,5}, {2,4,6}) | draft |
| D4 | draft Sec. 6.7 sign-flip print header | draft (0.4.2 text) |
| D5 | draft Secs. 2.7, 3.2, 6.5: `trim = "levels"` and the old irregular note | draft (0.4.2 text) |
| D6 | draft Secs. 2.7, 3.7, 6.7: default n_flip = 8 (now 6) | draft (0.4.2 text) |
| D7 | GTW Sec. 6.3 (and the draft): layout validity under a shared zeta_l | paper and draft; the package docs already hedge |
| D10 | CI end points are closures (documented, conservative) | none needed; any change is numeric and needs sign-off |

The four "known leads" in the task brief were all verified:
* Algorithm 1's mod: D1, and the example's inverse order: D2;
* draft Sec. 3.6: D3. `block_size` is K + 1, the block length in
  relabelled coordinates and the group order;
* rank instead of N - 2p: the code uses the numerical rank, confirmed
  against the reference on rank-deficient designs;
* the even-rep median: the reported value can leave the grid. The draft
  (Sec. 2.5) and `print()` both say so, so this is not a discrepancy.

#### Mutation check

| Mutation | testthat tests failing | base-R files failing | caught by (testthat, examples) |
|---|---:|---|---|
| (a) `<=` becomes `<` in Eq. (10) (`.ipt_eval()` and `.pval_matrix()`) | 3 | test-golden.R, test-irregular.R, test-lower-level.R, test-validation.R | test-draft-replication.R :: Sec. 6.5: the irregular example's numbers; test-pvalue-properties.R :: ties count toward the p-value (<=, not <); test-pvalue-properties.R :: D in col(X): p = 1 exactly, as Eq. (10) gives (and a warning) |
| (b) min_j a_j becomes a_1 (both places) | 15 | test-dyadic.R, test-golden.R, test-lower-level.R, test-panel.R, test-readme.R, test-signflip.R | test-designs.R :: three-way: independent groups per dimension, applied jointly; test-designs.R :: panel: p-values match Procedure 1 with period dummies in X; test-draft-replication.R :: Sec. 3.4: the first example prints as in the draft; +12 more |
| (c) project out col(X) only, not col([X, X_k]) | 19 | test-dyadic.R, test-equivariance.R, test-golden.R, test-lower-level.R, test-panel.R, test-readme.R, test-signflip.R | test-designs.R :: three-way: independent groups per dimension, applied jointly; test-designs.R :: panel: p-values match Procedure 1 with period dummies in X; test-designs.R :: irregular (Section 6.4): each repetition cuts cells to L0 at random; +16 more |
| (d) a fresh (pi, sigma) in every panel period | 2 | test-golden.R, test-panel.R, test-readme.R | test-designs.R :: panel: p-values match Procedure 1 with period dummies in X; test-draft-replication.R :: Sec. 6.2: the panel fit, with and without time effects |
| (e) K independent random permutations in place of Algorithm 1's cyclic group | 17 | test-dyadic.R, test-golden.R, test-lower-level.R, test-panel.R, test-readme.R, test-signflip.R | test-algorithm1.R :: closed under composition and inverses (Proposition 2); test-algorithm1.R :: psi_k is psi_1 composed with itself k times; test-algorithm1.R :: every non-identity element moves (K+1) floor(n/(K+1)) indices; +14 more |

Each mutation was applied alone to a throwaway worktree on branch `mutation-tmp`. The testthat suite ran from source, and the base-R suite against a private install of the mutant. The file was restored before the next mutation. The worktree and the branch were then deleted. Every mutation is caught by both suites.

Two notes on the pattern:

* Mutation (a) is caught only where ties are real (y = 0, a cell-constant d, D in col(X), and the irregular example, whose cell-constant d gives exact ties). The reference comparisons avoid near-ties on purpose.
* Mutation (d) is caught by the two p-value comparisons with the shared-(pi, sigma) reference, and by the golden and panel base files, but not by the structural test. That test cannot reach the front end's own group; see F-1.

#### Runtimes

About 30 s for the default testthat suite; 140 s for the slow Monte Carlo
tests alone; about 60 s for the base-R suite; 155 s for
`R CMD check --as-cran` (all tests included).

#### Decisions made without the maintainer (please review)

1. **A new branch, and a snapshot commit of your uncommitted 0.4.2 work.**
   The brief asked for small commits on a new branch. The 0.4.2 changes
   sat uncommitted in the same files, so they were committed unchanged as
   `1b62733` ("Snapshot: ...") on `test-suite-paper-fidelity`;
   `release/0.4.1` is untouched.
   * To undo while keeping everything uncommitted:
     `git switch test-suite-paper-fidelity && git reset --soft 4d166c9`.
   * Nothing was pushed.
2. **testthat is alongside, not instead.** The base-R suite was kept rather
   than converted: converting means deleting files this session did not
   create.
3. **gravity is not declared.** The Sec. 7 test holds the package name in a
   variable, so R CMD check does not demand it in Suggests (which would make
   every CI leg install it). If you prefer it declared, add
   `gravity` to Suggests and use the literal name.
4. **Sources.**
   * The JSS draft read was the 2026-09-18 build, outside the repository,
     located from the session notes. Its section numbering differs from
     the brief's by one after Sec. 6.4.
   * The RPT paper was read from page images, since its text layer is
     letter-shifted.
5. **Tools.** devtools is not installed. `testthat::test_local()`,
   `R CMD check --as-cran` and `roxygen2::roxygenise()` stood in for it.
   covr and roxygen2 7.3.2 came from an existing local tool library, so
   nothing was downloaded.

#### Follow-ups

See *Follow-ups (testing seams, not implemented)* above: F-1 to F-6.
