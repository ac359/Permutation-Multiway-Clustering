# `mwperm` test suite

The tests ship with the package: `tests/` is tracked in git, included in the
built tarball, and run by `R CMD check`. They are plain base-R `stopifnot()`
scripts — **`testthat` is not a dependency and is not used** — and they are
deliberately kept fast, because every one of them runs on every check.

## Running them

Install first, then run the top-level files:

```bash
R CMD INSTALL . && for f in tests/test-*.R; do Rscript "$f" || echo "FAILED: $f"; done
```

Each file prints `<file>: all assertions passed` when it finishes.

> **Install first — this has cost real time twice.** Many files reach package
> internals through `mwperm:::`, which resolves against the **installed**
> package, never against the working tree. Running a test after editing `R/`
> without reinstalling silently exercises the *old* code, or fails with a
> confusing "object not found". `devtools::load_all()` and `source()`-ing `R/`
> are fine for interactive poking at one function, not for trusting a result.

To run a single lower-level file on its own:

```bash
Rscript tests/lower-level-tests/test-permset.R
```

## Layout

```
tests/
  test-*.R              one file per user-facing entry point
  lower-level-tests/    the machinery underneath
  helpers/assertions.R  shared assertion helpers (sourced, never run as a test)
  golden/               seeded reference snapshots and their generator
```

`R CMD check` runs only the `.R` files at the **top level** of `tests/`, so the
lower-level suite is driven by `test-lower-level.R`, which sources every file in
`lower-level-tests/` in a fresh environment. A failure down there is a failure
here; nothing is skipped in a check.

### Top level — one file per method

Each file answers *"does this function do what its help page promises?"* and is
named after the function it tests.

| File | Covers |
|---|---|
| `test-main.R` | `mwperm()`, `mwperm_check()` — design detection, the printed diagnosis, forced-design validation, the assumption warnings, and **dispatch identity**: `mwperm()` must return exactly what the direct front-end call returns — including with `n_reps = NULL` (the front end's own default: 500 for irregular, 10 elsewhere; 0.4.2) and with `L0` forwarded to the incomplete-panel design |
| `test-formula.R` | `mwperm_formula()` — the `y ~ d \| x` interface; an identity test, since it only assembles arguments |
| `test-dyadic.R` | `mwperm_dyadic()` — the **two-way** (dyadic) design, condition InvA |
| `test-signflip.R` | `mwperm_dyadic_het()` and `build_flip_set()` — the same two-way design under **sign symmetry** instead of exchangeability (the heteroskedasticity-robust sign-flip test, 0.4.0): agreement with a corrected port of the author's `RPT_signflip()` given the same flip-group assignment (`helpers/signflip-reference.R`), the group axioms and the `2^(n_flip-1)` order, the `.apply_op()` signed-gather contract, size under heteroskedasticity, the `n_flip >= 6` resolution guard, seed hygiene, the opt-in dispatcher route, that the permutation designs are untouched, and (0.4.1) incomplete arrays: the `cells` kernel guard of `build_flip_set()` is bit-identical on a complete array and gives `2^(n_flip-1)` distinct elements on a diagonal-only one, and the fit on an incomplete array equals a from-scratch Procedure 1 with an explicit orthonormal complement. Named for the *group* rather than the function because `build_flip_set()` is tested here too |
| `test-threeway.R` | `mwperm_threeway()` — the three-way design, InvA in three dimensions |
| `test-panel.R` | `mwperm_panel()` — the panel design, condition InvB, including invariance to an arbitrary common time trend |
| `test-panel-missing.R` | `mwperm_panel_missing()` — incomplete panels; pins that on a complete array it builds permutations *identical* to `mwperm_panel()`. Section 8 (0.4.2) is `L0`, the same L0 periods in every cell: `L0 = n_t` equals `L0 = NULL` number for number, every gather vector of the front end's own builder preserves the period, every retained pair observes all of S, `L0 = 2` on a 4-period array equals the `L0 = NULL` fit on the data subset to those periods (these assertions moved here from the irregular file's retired `trim = "levels"` section) |
| `test-layout.R` | `mwperm_layout()` — replicated two-way layouts, within-cell exchangeability, `L0` balancing |
| `test-missing.R` | `mwperm_missing()` and `find_bicliques()` — incomplete arrays via fully observed blocks, one-sided permutation |
| `test-irregular.R` | `mwperm_irregular()` — Section 6.4 step (i) as printed (the only path since 0.4.2; `n_reps` defaults to 500 and is passed explicitly throughout): on a balanced array with `L0 = T` the gather vectors equal the panel construction's; every repetition keeps exactly `L0` observations per retained cell, redrawn from the rep seed (and reproduced by the rep-parallel path); every element holds the within-cell position fixed; the `N > 2p` check counts `cells_used * L0`; `trim` is refused by `mwperm()` with the replacement named; every fit's first note states that the trim is exact only for exchangeable replicates and names `mwperm_panel_missing()`; a second note fires when `d` varies inside a retained cell and not otherwise. Section 8 pins the engine's per-repetition `"rows"` hook against `mwperm_dyadic()` on a subset |
| `test-methods.R` | `print` / `summary` / `confint` / `coef` / `nobs` / `plot` / `mwperm_save` — the **output-label contract** (including the print header, which names the sign-flip group for `mwperm_dyadic_het()` fits) and every figure path |
| `test-validation.R` | the input-validation contract shared by all front ends: bad input fails early naming the argument, degenerate input gets the exact answer, and no confidence set is reported that the design cannot support |
| `test-equivariance.R` | properties every fit must have: row-order, nuisance (FWL), null-shift and scale/sign equivariance, and **parallel ≡ serial**, bit for bit (including the sign-flip, irregular and `L0` incomplete-panel builders) |
| `test-readme.R` | the transcripts `README.md` shows are the ones the package prints |
| `test-golden.R` | the seeded snapshot gate — 30 entries across all eight designs compared field by field against `golden/baseline.rds` |
| `test-lower-level.R` | the runner for `lower-level-tests/` |

Design-specific arguments are tested in that design's own file (`L0` and `rep`
in `test-layout.R`, `permute` and `min_block` in `test-missing.R`, `time_fe` in
`test-panel.R`); everything shared is in `test-validation.R`,
`test-methods.R` or `test-equivariance.R`.

`test-readme.R` is the one file that reaches outside the package: it locates
`README.md` by content, because `tests/` has a README of its own. Under
`R CMD check` the package README is not shipped alongside the tests, so that
half skips; the half that pins the printed lines always runs.

### `lower-level-tests/` — the machinery underneath

Each file answers *"is the mathematical object the method is built from
correct?"*. These reach into package internals and are keyed to the equations
of Guo, Toulis & Wang (2026) (GTW) and Wen, Wang & Wang (2025).

| File | Covers |
|---|---|
| `test-permset.R` | `build_perm_set()` — GTW Algorithm 1: a cyclic group of order `K+1`, identity first, closed under composition, RNG-hygienic |
| `test-obsperms.R` | the observation-level gather vectors for every design; the one-sided block subgroups; unobserved-cell detection; the `.cell_code` exactness guard |
| `test-bijection.R` | every gather vector is a bijection — a duplicated cell would silently make it many-to-one |
| `test-projection.R` | the FWL / partialling step (GTW Eq. 3) against an **independent** from-scratch projector that shares no code with `R/` |
| `test-pvalue.R` | the minorized randomization p-value (GTW Eq. 10): the grid it lives on, conservative ties, the `1/(K+1)` floor |
| `test-aggregate.R` | the cross-repetition rule (GTW Remark 1): one aggregation for the p-value *and* the confidence set |
| `test-exact-ci.R` | the closed-form confidence set: breakpoints, disconnected components, and what the reported end points mean |
| `test-invert-ci-grid.R` | the bracketing fallbacks: explicit-`grid` inversion, region budget guards, the disconnected-acceptance island guard |

## Where does a new test go?

* Testing what an **exported function** does for a user → the top-level file
  named after that function.
* Testing an **internal** — a permutation group, a projector, a p-value, an
  inversion routine → `lower-level-tests/`.
* Adding a **new design** → a new top-level `test-<design>.R`, plus a row in
  the table above. A new statistical feature also needs a null size check and a
  power sanity check; those live in `inst/replication/`, not here, because they
  are too slow for `R CMD check`.

## `helpers/`

`assertions.R` holds `msg_of()`, `expect_err()`, `expect_warn()`, `warns_of()`,
`msgs_of()`, `same_fit()`, `internal()` and `passed()`; `signflip-reference.R`
holds the corrected port of the author's `RPT_signflip()` and the
heteroskedastic gravity DGP that `test-signflip.R` checks the package
against. They live in a subdirectory so `R CMD check` does not mistake them
for tests, and every file loads what it needs with the same two-line idiom,
which resolves both from the package root and from `tests/`:

```r
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))
```

## `golden/`

`baseline.rds` is the authoritative record of every seeded number the package
produces; `baseline-0.4.0.rds` and `baseline-0.2.0.rds` are earlier releases'
snapshots, kept so the change lists in `NEWS.md` can be re-derived
(`make_baseline.R --check --against=baseline-0.4.0.rds`). Regenerate only after an intentional,
documented change:

```bash
R CMD INSTALL . && Rscript tests/golden/make_baseline.R
```

A drift in `pvalue`, `estimate`, `conf_int` or `K` is a defect until proven
otherwise — the paper's tables and the README's shown output are keyed to these
numbers.
