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
| `test-main.R` | `mwperm()`, `mwperm_check()` — design detection, the printed diagnosis, forced-design validation, the assumption warnings, and **dispatch identity**: `mwperm()` must return exactly what the direct front-end call returns |
| `test-formula.R` | `mwperm_formula()` — the `y ~ d \| x` interface; an identity test, since it only assembles arguments |
| `test-dyadic.R` | `mwperm_dyadic()` — the **two-way** (dyadic) design, condition InvA |
| `test-threeway.R` | `mwperm_threeway()` — the three-way design, InvA in three dimensions |
| `test-panel.R` | `mwperm_panel()` — the panel design, condition InvB, including invariance to an arbitrary common time trend |
| `test-layout.R` | `mwperm_layout()` — replicated two-way layouts, within-cell exchangeability, `L0` balancing |
| `test-missing.R` | `mwperm_missing()` and `find_bicliques()` — incomplete arrays via fully observed blocks, one-sided permutation |
| `test-irregular.R` | `mwperm_irregular()` — Section 6.4, irregular within-cell counts |
| `test-methods.R` | `print` / `summary` / `confint` / `coef` / `nobs` / `plot` / `mwperm_save` — the **output-label contract** and every figure path |
| `test-validation.R` | the input-validation contract shared by all front ends: bad input fails early naming the argument, degenerate input gets the exact answer, and no confidence set is reported that the design cannot support |
| `test-equivariance.R` | properties every fit must have: row-order, nuisance (FWL), null-shift and scale/sign equivariance, and **parallel ≡ serial**, bit for bit |
| `test-golden.R` | the seeded snapshot gate — 22 fits across all six designs compared field by field against `golden/baseline.rds` |
| `test-lower-level.R` | the runner for `lower-level-tests/` |

Design-specific arguments are tested in that design's own file (`L0` and `rep`
in `test-layout.R`, `permute` and `min_block` in `test-missing.R`, `time_fe` in
`test-panel.R`); everything shared is in `test-validation.R`,
`test-methods.R` or `test-equivariance.R`.

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

## `helpers/assertions.R`

`msg_of()`, `expect_err()`, `expect_warn()`, `warns_of()`, `msgs_of()`,
`same_fit()`, `internal()` and `passed()`. It lives in a subdirectory so
`R CMD check` does not mistake it for a test, and every file loads it with the
same two lines, which resolve both from the package root and from `tests/`:

```r
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))
```

## `golden/`

`baseline.rds` is the authoritative record of every seeded number the package
produces; `baseline-0.2.0.rds` is the previous release's, kept so the change
list in `NEWS.md` can be re-derived. Regenerate only after an intentional,
documented change:

```bash
R CMD INSTALL . && Rscript tests/golden/make_baseline.R
```

A drift in `pvalue`, `estimate`, `conf_int` or `K` is a defect until proven
otherwise — the paper's tables and the README's shown output are keyed to these
numbers.
