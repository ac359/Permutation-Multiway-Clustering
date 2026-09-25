# `mwperm` test suite

Everything under `tests/` is tracked in git, ships in the package tarball and
runs under `R CMD check`, so every file here is kept fast.

## Two suites, two questions

| | Base-R suite | testthat suite |
|---|---|---|
| **Question it answers** | Does each exported function keep its contract, and do the seeded numbers stay put? | Does the code compute what the paper defines? |
| **Where** | `tests/test-*.R` and `tests/lower-level-tests/` | `tests/testthat/` |
| **Organised by** | exported function (one file each), then the machinery underneath | paper object (Algorithm 1, Procedure 1, Eq. 10, ...) |
| **Checked against** | the package's own contract and the golden seeded snapshot | a naive Procedure 1 written from the paper (`helper-reference.R`), the paper's algebra, the JSS draft's printed output, Monte Carlo |
| **Style** | plain `stopifnot()` scripts | testthat 3rd edition |
| **Tests which copy of the package** | the **installed** one: reinstall first | the source tree (`test_local()`) |
| **Time** | about 60 s | about 30 s, plus 2.5 min for the slow Monte Carlo tests |

The base-R suite came first and is the package's regression contract. The
testthat suite was added in 0.4.2 to show paper fidelity; it sits alongside
the base-R suite rather than replacing it, so the contract was not rewritten.
Where the two touch the same object they check it from different sides; the
[crosswalk](#where-each-topic-is-tested) lists both. `testthat (>= 3.0.0)` is in
`Suggests` for the second suite only.

## Running

From the package root:

```bash
# base-R suite: install first, then run the top-level files
R CMD INSTALL . && for f in tests/test-*.R; do Rscript "$f" || echo "FAILED: $f"; done

# testthat suite (loads the package from source; no install needed)
NOT_CRAN=true Rscript -e 'testthat::test_local()'

# ... plus the slow Monte Carlo checks
NOT_CRAN=true MWPERM_SLOW_TESTS=true Rscript -e 'testthat::test_local()'

# golden identity gate on its own
R CMD INSTALL . && Rscript tests/golden/make_baseline.R --check

# everything, as CRAN and CI run it
R CMD build . && R CMD check --as-cran mwperm_*.tar.gz
```

Each base-R file prints `<file>: all assertions passed`. To run one
lower-level file: `Rscript tests/lower-level-tests/test-permset.R`. To run
one testthat file: `testthat::test_local(filter = "algorithm1")`.

> **Install before trusting a base-R result.** These files reach package
> internals through `mwperm:::`, which resolves against the **installed**
> package, never the working tree. After editing `R/` without reinstalling,
> they silently test the *old* code, or fail with a confusing "object not
> found". This has cost real time twice. `devtools::load_all()` and
> `source()`-ing `R/` are fine for poking at one function, not for trusting a
> result.

Two environment variables widen the testthat suite:

* `MWPERM_SLOW_TESTS=true` runs `test-montecarlo.R` (Theorem 1 checked at
  every atom of the p-value distribution). It never runs on CRAN.
* `MWPERM_SHOW_DISCREPANCIES=true` runs the tests skipped with
  `Discrepancy D<n>`. Each asserts what the paper or the JSS draft says where
  the code deliberately does something else; setting the variable shows each
  one failing. The entries are explained in the discrepancy log of
  `TESTING_PLAN.md` at the repository root.

## Where does a new test go?

| You are testing... | Put it in |
|---|---|
| what an **exported function** promises: its seeded numbers, arguments, errors, notes | the base-R file named after that function (`test-dyadic.R` for `mwperm_dyadic()`, ...) |
| a contract **shared by every front end**: scalar-argument validation, output labels, equivariance, parallel identity | base-R `test-validation.R`, `test-methods.R` or `test-equivariance.R` |
| an **internal building block** (a group, the projector, the p-value, an inversion routine) against its own invariants | base-R `lower-level-tests/`, one file per object |
| that the code **matches a paper object** (against the naive reference, the paper's algebra, or a printed draft output) | the testthat file named after that paper object |
| a place where the code **deliberately differs** from the paper or the draft | a testthat test asserting the paper's version, skipped with `skip_discrepancy("D<n>", ...)`, plus an entry in the discrepancy log of `TESTING_PLAN.md` |
| **size or power** by simulation | quick checks: `testthat/test-montecarlo.R` (gated by `MWPERM_SLOW_TESTS`); publication-scale tables: `inst/replication/`, which no check runs |

**A new design** needs all of these:
- a base-R `test-<design>.R` and a row in the table below;
- an entry in `golden/make_baseline.R`, then a regenerated baseline;
- a reference check in `testthat/test-designs.R`, comparing its p-values with
  the naive Procedure 1 on the group rebuilt from the seed;
- a size check in `testthat/test-montecarlo.R`;
- a row in the traceability table of `TESTING_PLAN.md`.

## Where each topic is tested

`ll/` is `lower-level-tests/`.

| Topic | Base-R suite | testthat suite |
|---|---|---|
| Algorithm 1, `build_perm_set()` | `ll/test-permset.R` | `test-algorithm1.R` |
| each design's group, lifted to observations | `ll/test-obsperms.R`, `ll/test-bijection.R`, each design's file | `test-designs.R` |
| Procedure 1: projection and statistic | `ll/test-projection.R` | `test-procedure1-reference.R` |
| the p-value, Eq. 10 | `ll/test-pvalue.R` | `test-pvalue-properties.R` |
| aggregation over repetitions (Remark 1) | `ll/test-aggregate.R` | `test-pvalue-properties.R` |
| confidence sets by inversion | `ll/test-exact-ci.R`, `ll/test-invert-ci-grid.R` | `test-confidence-sets.R` |
| invariance and equivariance | `test-equivariance.R` | `test-invariances.R` |
| seeds, parallel equals serial, the caller's RNG | `test-equivariance.R`, each design's file | `test-reproducibility.R` |
| `mwperm()` dispatch and `mwperm_check()` | `test-main.R` | `test-dispatch.R` |
| S3 methods and output labels | `test-methods.R` | `test-s3-methods.R` |
| input validation | `test-validation.R` (scalar arguments) | `test-data-validation.R` (data, interfaces) |
| `find_bicliques()` | `test-missing.R` | `test-bicliques.R` |
| the sign-flip test | `test-signflip.R` | `test-dispatch.R`, `test-reproducibility.R`, `test-pvalue-properties.R`, `test-draft-replication.R` |
| seeded numbers people see | `test-golden.R`, `test-readme.R` | `test-draft-replication.R` |
| validity by simulation | none (too slow; see `inst/replication/`) | `test-montecarlo.R` |

## Layout

```
tests/
  README.md                 this file
  test-*.R                  base-R suite: one file per exported function,
                              plus the shared contracts and the golden gate
  test-lower-level.R        runs every file in lower-level-tests/
  lower-level-tests/        base-R suite: the machinery underneath
  helpers/                  sourced by base-R files, never run as tests
  golden/                   the seeded snapshots and their generator
  testthat.R                runner that R CMD check uses for testthat/
  testthat/
    helper-reference.R      the naive Procedure 1 and shared generators
    test-*.R                testthat suite: one file per paper object
    _snaps/                 expect_snapshot() records
```

`R CMD check` runs only the `.R` files at the **top level** of `tests/`.
That is why the lower-level files are driven by `test-lower-level.R`, which
sources each in a fresh environment and lets any failure propagate, and why
`helpers/` is a subdirectory: it keeps them from being mistaken for tests.

## Base-R suite: top level, one file per exported function

Each file answers *"does this function do what its help page promises?"*,
and opens with a header saying what it pins and why.

| File | Covers |
|---|---|
| `test-main.R` | `mwperm()`, `mwperm_check()`: design detection, the printed diagnosis, forced designs, the assumption warnings, and **dispatch identity** (`mwperm()` returns exactly what the direct front-end call returns) |
| `test-formula.R` | `mwperm_formula()`: an identity test against the data interface, since it only assembles arguments |
| `test-dyadic.R` | `mwperm_dyadic()`: the two-way design, condition InvA |
| `test-threeway.R` | `mwperm_threeway()`: InvA in three dimensions |
| `test-panel.R` | `mwperm_panel()`: condition InvB, including invariance to an arbitrary common time trend; `time_fe` |
| `test-panel-missing.R` | `mwperm_panel_missing()`: incomplete panels; on a complete array its permutations are *identical* to `mwperm_panel()`'s; `L0` |
| `test-layout.R` | `mwperm_layout()`: within-cell permutations, `rep`, `L0` balancing |
| `test-missing.R` | `mwperm_missing()` and `find_bicliques()`: fully observed blocks, `permute`, `min_block` |
| `test-irregular.R` | `mwperm_irregular()`: Section 6.4 step (i), the per-repetition subsample, its notes, and the engine's `"rows"` hook |
| `test-signflip.R` | `mwperm_dyadic_het()` and `build_flip_set()`: agreement with a corrected port of the authors' `RPT_signflip()` (`helpers/signflip-reference.R`), the group's order `2^(n_flip-1)`, the `.apply_op()` contract, size under heteroskedasticity, incomplete arrays |
| `test-methods.R` | `print`, `summary`, `confint`, `coef`, `nobs`, `plot`, `mwperm_save`: the **output-label contract** and every figure path |
| `test-validation.R` | scalar-argument validation shared by all front ends: bad input fails early and names the argument; degenerate input gets the exact answer; no confidence set the design cannot support |
| `test-equivariance.R` | row order, nuisance (FWL), null shift, scale and sign equivariance, and **parallel equals serial**, bit for bit |
| `test-readme.R` | the transcripts `README.md` shows are what the package prints |
| `test-golden.R` | the seeded snapshot gate: 30 entries across all eight designs, compared field by field with `golden/baseline.rds` |
| `test-lower-level.R` | the runner for `lower-level-tests/` |

`test-readme.R` is the one file that reaches outside the package. It finds
`README.md` by content, because `tests/` has a README of its own. Under
`R CMD check` the package README is not next to the tests, so that half
skips; the half that pins the printed lines always runs.

## Base-R suite: `lower-level-tests/`, the machinery

Each file answers *"is the mathematical object the method is built from
correct?"*. They reach into package internals and are keyed to the equations
of Guo, Toulis & Wang (2026) (GTW) and Wen, Wang & Wang (2025).

| File | Covers |
|---|---|
| `test-permset.R` | `build_perm_set()`, GTW Algorithm 1: a cyclic group of order `K+1`, identity first, closed under composition, RNG-hygienic |
| `test-obsperms.R` | the observation-level gather vectors for every design; the one-sided block subgroups; unobserved-cell detection; the `.cell_code()` exactness guard |
| `test-bijection.R` | every gather vector is a bijection; a duplicated cell would silently make it many-to-one |
| `test-projection.R` | the partialling step (GTW Eq. 3) against an independent projector that shares no code with `R/` |
| `test-pvalue.R` | the minorized p-value (GTW Eq. 10): its grid, conservative ties, the `1/(K+1)` floor |
| `test-aggregate.R` | the cross-repetition rule (GTW Remark 1): one aggregation for the p-value *and* the confidence set |
| `test-exact-ci.R` | the closed-form confidence set: breakpoints, disconnected components, what the reported end points mean |
| `test-invert-ci-grid.R` | the fallbacks: explicit-`grid` inversion, region budget guards, the island guard |

## testthat suite: `tests/testthat/`, one file per paper object

| File | Covers |
|---|---|
| `helper-reference.R` | shared code, loaded by testthat before every file, never a test: a naive Procedure 1 (QR residuals with numerical rank, then `a_k`, `b_k` and Eq. 10); each front end's group rebuilt from its seed by explicit key matching; data generators; Monte Carlo helpers; `skip_discrepancy()` |
| `test-procedure1-reference.R` | the package's `a_k`, `b_k` (to 1e-8) and every `pvalues_rep` entry against the naive Procedure 1, including d > 1, a non-zero null, shuffled rows and rank-deficient `[X \| X_k]` |
| `test-algorithm1.R` | `build_perm_set()`: bijections, closure and inverses (Proposition 2), `psi_k = psi_1^k`, orbit sizes, the printed formula; the draft's Sec. 3.6 example |
| `test-pvalue-properties.R` | the p-value grid, the median over repetitions (odd and even), `median2`, ties counted (`<=`), a `D` the group cannot move, `D` in `col(X)` |
| `test-invariances.R` | `y + X gamma`, rescaling `y` or `D`, testing `b` versus `y - D b`, a common trend under `time_fe`; CI equivariance |
| `test-confidence-sets.R` | duality just inside and outside every end point, the resolution guard, bisection and grid against the exact set, the joint region |
| `test-designs.R` | each design's group structure (three-way, panel, layout and `L0`, missing, incomplete panel, irregular) and its p-values against the reference |
| `test-bicliques.R` | fully observed, disjoint, `min_block`; the exact search against brute force; the diagonal-deleted 40 x 40 array |
| `test-dispatch.R` | every row of the draft's Table 1; `mwperm()` identical to its worker; `mwperm_check()` snapshots |
| `test-s3-methods.R` | the draft's Sec. 3.5 slots, `summary()`, `confint()`, `coef()`, `nobs()`, `print()` snapshots, `plot()` |
| `test-data-validation.R` | lengths, missing and non-finite values, `N > 2p`, duplicate cells, the vector, column-name and formula interfaces, cluster-id types |
| `test-reproducibility.R` | same seed, `seed + r - 1`, `n_cores = 2` identical to serial, the global RNG state |
| `test-montecarlo.R` | slow: the empirical CDF at every atom `j/(K+1)` under the null (t_3, Cauchy, Moulton, MCAR, trending panel), the draft's negative control, and power |
| `test-draft-replication.R` | every seeded output the JSS draft prints, at printed precision; OLS against `lm()`; Sec. 7 when `gravity` is installed |

`_snaps/` holds the `expect_snapshot()` records for `test-dispatch.R` and
`test-s3-methods.R`. They are skipped on CRAN. After an intended change to
printed output, update them with `testthat::snapshot_accept()`.

## `helpers/`

`assertions.R` holds `msg_of()`, `expect_err()`, `expect_warn()`,
`warns_of()`, `msgs_of()`, `same_fit()`, `internal()` and `passed()`.
`signflip-reference.R` holds the corrected port of the authors'
`RPT_signflip()` and the heteroskedastic gravity DGP that `test-signflip.R`
checks against. Every base-R file loads what it needs with the same two
lines, which work from the package root and from `tests/`:

```r
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))
```

## `golden/`

`baseline.rds` is the authoritative record of every seeded number the
package produces, and `make_baseline.R` generates and checks it.
`baseline-0.4.1.rds`, `baseline-0.4.0.rds` and `baseline-0.2.0.rds` are
earlier releases' snapshots, kept so the change lists in `NEWS.md` can be
re-derived:

```bash
Rscript tests/golden/make_baseline.R --check --against=baseline-0.4.1.rds
```

Regenerate `baseline.rds` only after an intentional change that `NEWS.md`
documents:

```bash
R CMD INSTALL . && Rscript tests/golden/make_baseline.R
```

`--check` names the fields that differ. A difference only in `note` is a
wording change. A drift in `pvalue`, `estimate`, `conf_int` or `K` is a
defect until proven otherwise: the paper's tables and the README's shown
output are keyed to these numbers.
