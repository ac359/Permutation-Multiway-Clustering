# Paper → code map

Every object the method paper defines, and where `mwperm` implements it.

**Method paper.** Guo, W., Toulis, P. and Wang, Y. (2026). *Permutation
Inference under Multi-way Clustering and Missing Data*. arXiv:2601.08610
[stat.ME]. Authors, per the title page: Wenxuan Guo (Booth, Chicago), Panos
Toulis (Booth, Chicago), Yuhao Wang (IIIS, Tsinghua and Shanghai Qi Zhi).
Section, equation, procedure, algorithm and theorem numbers below are that
paper's.

**Companion source.** Wen, K., Wang, T. and Wang, Y. (2025). Residual
permutation test for regression coefficient testing. *Annals of Statistics*
**53**(2), 724–748. `mwperm`'s residual projection and the block-cyclic
permutation construction are the dyadic adaptations of that paper's Algorithm 1
and residual permutation test.

Read the "Deliberate deviation" column carefully: where it is not "—", the code
and the paper's printed text differ, and the reason is given. Every such
deviation is either (a) a case the paper's text does not cover and the code must
resolve, or (b) a place where the printed formula is internally inconsistent and
the mathematically coherent reading was implemented instead. None of them
changes the object the theorems are about.

---

## Model, assumptions, and the null

| Paper object | Where | Deliberate deviation |
|---|---|---|
| Model (1)/(7): `y = Xγ + Dβ + ε` stacked over cells | assembled per design in the front ends (`R/dyadic.R`, `R/panel.R`, `R/threeway.R`, `R/layout.R`, `R/irregular.R`, `R/missing.R`); the intercept is added by `.make_X()` in `R/engine.R` | — |
| `H0: β = b` (the paper states `β = 0`) | `beta_null` argument, applied in `.ipt_eval()` (`R/core.R`) | Generalised from `β = 0` to any fixed `b`. The shift is applied to **both** terms — `a_k` uses `y − Db`, `b_k` uses `y_k − D_k b` — which is the exact analogue and is what makes test inversion coherent. Applying it to only one term would silently break the duality. |
| **Assumption 1** — conditional separate (double) exchangeability of the errors, `(ε_ij) =d (ε_π(i)σ(j)) | X, D` | the validity premise of `mwperm_dyadic()` (`R/dyadic.R`) and, in its three-way form (InvA, §6.1), `mwperm_threeway()` (`R/threeway.R`) | — |
| **Assumption 4** — the missingness mask is independent of the errors given the covariates | the validity premise of `mwperm_missing()` (`R/missing.R`) and of `mwperm_irregular()` (`R/irregular.R`), where the "mask" is `M_ij = 1{ℓ_ij ≥ L0}`. Stated in both roxygen blocks. | — |
| Assumptions 2, 3, 5 (moment and random-graph conditions) | not implemented — they are conditions for the *power* theorems, not the test | — |

## Procedure 1 — the invariant permutation test

| Paper object | Where | Deliberate deviation |
|---|---|---|
| **Procedure 1, step 1** — orthonormal `V_k` with `V_k'X = 0` and `V_k'X_{π_k,σ_k} = 0`; compute `a_k = ‖D'V_kV_k'y‖`, `b_k = ‖D'V_kV_k'y_{π_k,σ_k}‖` | `.ipt_prepare()` (`R/core.R`) caches `u = Dr'y`, `v = Dr'y_k`, `M = Dr'Dr`, `W = Dr'D_k` where `Dr = V_kV_k'D`; `.ipt_eval()` reassembles `a_k(b)`, `b_k(b)` | **Three.** (i) `V_k` is never materialised: `V_kV_k' = I − P_{[X|X_k]}` is symmetric idempotent, so `D'V_kV_k'z = (V_kV_k'D)'z` and only `D` is residualised — same numbers, one projection instead of two. (ii) The projector goes through an **eigendecomposition pseudo-inverse of the 2p×2p Gram matrix**, not an explicit `qr.Q`. Both diagonal Gram blocks equal `X'X` because a permutation only reorders rows; the rank cut is `1e-14` because eigenvalues of a Gram matrix are *squared* singular values (`.lm.fit`'s `1e-7`, squared). (iii) **The projection dimension is `N − rank([X | X_k])`, not the paper's `N − 2p`.** The stack is rank-deficient *by construction* — a permutation maps the intercept to itself, and in `mwperm_panel(time_fe = TRUE)` maps the period dummies among themselves — so `N − 2p` is not even well defined here. The two orthogonality conditions are what the proof uses, and they are satisfied exactly. This is also a deliberate difference from the paper authors' reference code, whose `get_orthogonal_matrix` takes the last `n − ncol(X_aug)` columns of a full `qr.Q` and so assumes full column rank. |
| **Equation (10)** — `pval = (1 + Σ_k 1{min_{1≤j≤K} a_j ≤ b_k}) / (K+1)` | `.ipt_eval()` (`R/core.R`): `(1 + sum(b >= amin)) / Kp1` | — . The minorisation `min_j a_j` runs over the **non-identity** permutations `j = 1..K` only, and the comparison keeps the paper's tie direction (`b_k ≥ min_j a_j` counts). Both are load-bearing: reversing the tie would make the test anti-conservative at the boundary. The reject rule is `p ≤ alpha`, and p lives on the grid `{1,…,K+1}/(K+1)`. |
| **Procedure 1, step 3** — `CI = {x : pval(x) > α}` | `.exact_ci_set()` and `.invert_ci()` (`R/engine.R`) for `d = 1`; `.invert_region()` for `d > 1` | For `d = 1` the set is computed **exactly**, in closed form, rather than approximated. `a_j(b)` and `b_k(b)` are piecewise linear, so `pval` is a step function whose jumps solve `|v_k − W_k b| = |u_j − M_j b|`, i.e. `b = (u_j − v_k)/(M_j − W_k)` and `b = (u_j + v_k)/(M_j + W_k)`; the code evaluates the aggregated p-value at every root and inside every cell between roots. No bracketing assumption, no bisection tolerance, and a disconnected set is reported as such (`$conf_set`) instead of being replaced by its hull. `$conf_int` remains the hull for backward compatibility. Above a candidate budget (≈ `2K²·n_reps`) it falls back to bracket-and-bisect and says so in `$note`. |
| **Theorem 1** — `P(pval ≤ α | X, D) ≤ α` for all α, all n | the guarantee the package sells; exercised by the size and coverage simulations in `inst/replication/` | — . Theorem 1 covers **one** random group, i.e. `n_reps = 1`. See Remark 1 below for what `n_reps > 1` does and does not inherit. |
| **Remark 1** — run the procedure several times and take the median p-value | `.agg_pvals()` (`R/engine.R`), applied in exactly two places: the reported p-value and the confidence set | The median of dependent randomised p-values is a de-randomisation heuristic, not a theorem: it is endorsed by Remark 1 but is not itself a level-α p-value. `aggregate = "median2"` returns `min(1, 2 × median)`, which is (Rüschendorf 1982; Vovk and Wang 2020) and inverts the same rule. Default unchanged. **The rule is applied once and used for both the p-value and the set**, so the test and the interval cannot disagree — before 0.3.0 they could, because the default interval took the median of per-rep *end points* (not an inversion of anything) while the explicit-`grid` path took a union across reps. |
| **Theorem 2** — power, for **fixed** `K` with `n` divisible by `K+1` | not implemented (a theorem, not a routine); its regime is noted in `?mwperm_dyadic` | The **default** `K = min(n_row, n_col) − 1` is deliberately *not* Theorem 2's regime: it makes `K+1` equal to the cluster count, so Algorithm 1 produces a single block spanning all clusters with no fixed tail. That maximises p-value resolution (`1/(K+1)`), which is what a user needs to reject at all, and validity (Theorem 1) holds for any group Algorithm 1 builds. Set `K` explicitly to work in Theorem 2's regime. |
| Theorem 3 (power without random effects) | not implemented | — |

## Algorithm 1 — permutation-set construction

| Paper object | Where | Deliberate deviation |
|---|---|---|
| **Algorithm 1** — random relabelling `π`, consecutive blocks of size `K+1`, tail held fixed, `ψ_k = π⁻¹ ∘ ψ̃_k ∘ π` | `build_perm_set()` (`R/perm_set.R`); RNG hygiene via `.save_seed()`/`.restore_seed()` | **The printed formula for `ψ̃_k` is not implemented literally, because read literally it is not a permutation.** The paper gives `ψ̃_k(i) = i + k` if `i mod (K+1) ≤ K+1−k`, else `i − (K+1−k)`. At `i ≡ 0 (mod K+1)` — the last element of a block — the first branch applies (`0 ≤ K+1−k`) and maps `i` to `i + k`, outside its block; two elements then collide and the map is not one-to-one. The paper's own worked example describes the intended object unambiguously and in the *opposite* shift direction: `S₁^k := (K+2−k, …, K+1, 1, 2, …, K+1−k)`. The code implements a clean cyclic shift of each block, `new_pos = (pos + k) mod (K+1)`, which is a permutation, generates the same cyclic group of order `K+1`, and is closed under composition — the property the proof (Proposition 2) actually uses. Do **not** "fix" the code to the printed formula. |
| Group closure `π_k = π_r ∘ π_s` (Proposition 2) | guaranteed by construction — every element is a power of one generator; asserted in `tests/test-permset.R` | — |
| Applying Algorithm 1 twice (rows, columns) for dyadic regression | `mwperm_dyadic()` (`R/dyadic.R`), combined into observation gather vectors by `.build_obs_perms()` (`R/core.R`) | — |

## Procedure 2 and Algorithm 2 — missing data

| Paper object | Where | Deliberate deviation |
|---|---|---|
| **Procedure 2** — disjoint fully observed blocks; per-block row/column groups of **common order `K+1`**; concatenate into one group; apply Procedure 1 to the stacked data | `mwperm_missing()` (`R/missing.R`), with `.build_obs_perms_blocks()` building the block-diagonal gather vectors | One `K` is shared across all blocks and cells, per step 2 — enforced by `.default_K()` against the smallest permuted block side. Cells outside the selected blocks are discarded and the loss is reported in `$note`/`$cells_used`. `permute = "rows"`/`"cols"` adds a one-dimensional subgroup (the other margin held fixed) that the paper does not discuss; it is a subgroup of the same construction, so validity is unaffected, and it exists because a tall single-column block is fully usable when only rows are permuted. |
| **Theorem 4** — validity of Procedure 2 given Assumptions 1 and 4 | the guarantee for `mwperm_missing()` and `mwperm_irregular()` | — |
| **Algorithm 2** — iteratively solve the maximum-edge biclique problem and peel | `find_bicliques()` (`R/missing.R`), with `.grow_biclique()` (greedy, the default) and `.max_biclique_exact()` (branch and bound) | The paper proposes Bimax for the approximate solve; the package uses a seed-and-intersect greedy heuristic and an optional exact branch-and-bound with a node budget. The paper states explicitly that an approximate solve keeps Procedure 2 valid — a sub-maximal block costs power, never validity — so any heuristic is admissible here. Blocks are disjoint in **both** margins (`mwperm_missing()` relies on this; overlapping blocks must never be pooled). Under `method = "greedy"`, a peel that returns a block below `min_block` no longer ends the search immediately: it retries `retry_peels` times with the seed window slid down the degree order, because a small block from the highest-degree seeds is not evidence that no conforming block remains elsewhere. Under `method = "exact"` a sub-floor block *is* proof, and the search stops as before. |
| Remark 2 (graph-theoretic construction of `F_M`) | `find_bicliques()`, exported so the block structure can be inspected before fitting | — |
| Proposition 1, **Theorem 5** (largest biclique size; power under Assumption 5) | not implemented | — |

## Section 6 — extensions

| Paper object | Where | Deliberate deviation |
|---|---|---|
| **§6.1 Random effects / three-way (InvA)** — apply Algorithm 1 three times, then Procedure 1 | `mwperm_threeway()` (`R/threeway.R`) | — |
| **§6.2 Panel models (InvB)** — apply Algorithm 1 twice, over `[m]` and `[n]`; the **same** `(π, σ)` in every period, time held fixed | `mwperm_panel()` (`R/panel.R`); the time dimension gets a `NULL` group in `.build_obs_perms()`, which is exactly "held fixed" | `time_fe = TRUE` (the default) adds period dummies to the nuisance design. They are invariant to the within-period permutation, so validity is untouched; they de-bias the estimate. This is an addition to the printed procedure, not a change to it. |
| **§6.3 Two-way layouts** — apply Algorithm 1 on `[ℓ_ij]` **within each cell**, concatenate as in Procedure 2 | `mwperm_layout()` (`R/layout.R`), with `.build_obs_perms_layout()` | `L0` is available here as a **power lever**, but the `L0` threshold itself is from **§6.4**, not §6.3. `mwperm_layout(L0 = )` balances the array (drop thin cells, uniformly downsample the rest) and then runs the §6.3 *within-cell* test. That hybrid is valid — a uniform random subset of exchangeable replicates is exchangeable — and raising `L0` raises the attainable `K`, but it is **not** the §6.4 procedure and does not cover what §6.4 exists for. Documented as such in `?mwperm_layout` and in `README.md`. |
| **§6.4 Irregular designs** — mask `M_ij = 1{ℓ_ij ≥ L0}`; Algorithm 2 on `M`; randomly drop `ℓ_ij − L0` per retained cell; Procedure 2 on what remains | `mwperm_irregular()` (`R/irregular.R`), using `find_bicliques()` for the mask, `.downsample_to_L0()` for the deletion, and `.build_obs_perms_blocks(slot = )` for the permutation | Implemented as written. The within-cell slot is **held fixed** (the same device §6.2 uses for time): cell `(i,j)` slot `l` maps to cell `(π(i), σ(j))` slot `l`. One departure from §B of the appendix: the random deletion is drawn **once**, from `seed`, before the `n_reps` loop, so `n_reps` averages over the permutation draw but not over the subsample. §B re-draws both. To average over both, run several seeds and take the median across runs — documented in `?mwperm_irregular`, "Reproducibility". |
| §6.4's recommendation to tune `L0` by grid search | not automated | Choosing `L0` from the outcome would invalidate the guarantee. `?mwperm_irregular` says to choose it from the cell sizes, which are ancillary under Assumption 4, and the fit reports how many cells and observations survived. |
| §7 simulations, §8 gravity application, §B trust-level application | `inst/replication/` (seeded Monte-Carlo scripts and reference outputs) | The gravity intervals differ from the authors' captured output because *their* code omits the intercept from the nuisance projection. `mwperm` follows the paper; see the `Deliberate deviation` cell for Procedure 1 step 1 above, which records the same difference on the projection itself. |

## Objects with no paper counterpart

These exist for users and add no statistical content.

| Object | Where |
|---|---|
| `mwperm()` — design detection and dispatch | `R/unified.R` |
| `mwperm_check()` — the diagnosis without computing | `R/unified.R` |
| `mwperm_formula()` — `y ~ d | x` interface | `R/formula.R` |
| `print` / `summary` / `confint` / `coef` / `nobs` | `R/methods.R` |
| `plot` / `mwperm_save()` | `R/plot.R` |
| `trade_dyadic`, `trade_panel` (synthetic) | `data/`, generated by `data-raw/make_data.R` |

## Further paper-vs-code observations (recorded, not acted on)

Found while writing this map. None changes a computed statistic, so per the
working rules the code is left alone.

1. **`n_reps` default is 10, the paper's exposition is single-run.** Every front
   end defaults to `n_reps = 10`, so out of the box the package reports a
   *median* p-value — Remark 1's stabilisation — rather than the single-group
   p-value Theorem 1 covers. This is documented now in `?mwperm_dyadic` and
   `?confint.mwperm`, and `aggregate = "median2"` restores a guarantee that
   holds as stated. Confidence: high. Not a defect, but a paper describing
   `mwperm` must not say "the reported p-value is the Theorem 1 p-value" without
   qualifying it.
2. **`K` is capped at 199.** The paper imposes no cap. The cap bounds work and
   costs only resolution beyond `1/200`. Confidence: high that it is harmless.
3. **The layout test's cross-cell replicate effect.** §6.3's stated error
   structure is `ε_ijl = η_ij + ζ_l + u_ijl` with `ζ_l` *shared across cells*,
   but the procedure draws an **independent** permutation per cell, which
   changes the cross-cell alignment of `ζ`. The within-cell invariance argument
   covers `ε_ijl = η_ij + u_ijl`; it does not obviously cover a shared `ζ_l`.
   Simulation found no measurable size effect, and `?mwperm_layout` states the
   condition it actually needs. Confidence: medium that the paper's example is
   loosely stated rather than the code being wrong. Worth a line in the paper.
4. **Sub-seed collisions at scale.** Rep `r` uses `seed + r − 1`; within a rep,
   slot `j` used `rep_seed × 1000 + j`, which collides once `j ≥ 1000` (layouts
   with ≥ 1000 cells, block designs with ≥ 250 blocks). Each rep stayed exactly
   valid — the seed only picks the relabelling — but two reps shared a
   relabelling, so the median was averaging fewer effective draws than it
   thought. Fixed in 0.3.0 by widening the stride *only* for designs that would
   collide. Confidence: high.
