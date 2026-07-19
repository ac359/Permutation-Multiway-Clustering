# mwperm replication material

Self-contained Monte-Carlo scripts that reproduce the headline claims behind
`mwperm`: that the invariant permutation test (IPT) is **valid in finite
samples** across every design, and that its size control does not come at the
cost of power. Everything here uses only base R and the installed `mwperm`
package — no external data, no extra dependencies.

Locate the installed copy with

```r
system.file("replication", package = "mwperm")
```

and run the scripts from a writable copy of that directory (they read
`mc_lib.R` by a relative path and write to `./out` and `./cache`).

## What each script shows

| Script | Claim reproduced |
|---|---|
| `01_size.R` | Type-I error under the null is at/below nominal for the dyadic, panel, and three-way designs, and the whole p-value distribution is super-uniform — while the classical OLS *t*-test that ignores the clustering over-rejects grossly (the motivating failure). |
| `02_power.R` | Rejection rises monotonically toward 1 as the true coefficient moves away from the null (the test is not conservative-because-dead). |
| `03_ci_coverage.R` | The inverted-test confidence interval covers the truth at least at its nominal level (coverage is inherited from test validity). |
| `04_negative_control.R` | On a *trending* panel the three-way test (which permutes time) over-rejects badly while the panel test (which holds time fixed) stays valid — the reason `mwperm()` defaults to panel for 3-index arrays. |
| `05_permute.R` | The one-sided `permute = "rows"` option is valid and powerful on a design (a single fully observed column) the two-sided test cannot handle at all. |

## Running

With `mwperm` installed, from within a writable copy of this directory:

```sh
Rscript make.R                 # runs 01–05 in order, ~10–15 min total
MC_N=10000 Rscript make.R      # tighter Monte-Carlo error (much longer)
```

or run any script on its own, e.g. `Rscript 01_size.R`. Each script prints a
table and a pass/fail verdict, and writes it to `./out/<name>.txt` plus an
`.rds` summary. `make.R` also records `out/sessionInfo.txt`.

**Caching and resumability.** Every script is batched and checkpointed under
`./cache` (override with the `MWPERM_REPL_CACHE` environment variable). A
completed batch is never recomputed, so an interrupted run resumes where it
stopped — just re-run the same command. A cache directory keyed to a cell also
stores the parameters that produced it and refuses to be reused under
different parameters.

## Reference outputs

`expected/` holds the tables this suite produces at the default sim counts
(2000 per size/permute cell, 1500 per power cell, 600 per coverage cell, 1000
per negative-control cell), together with the `sessionInfo()` behind them.
Compare your
`out/` against `expected/` to confirm a faithful reproduction. Small
Monte-Carlo drift is expected — the numbers are random — but the qualitative
verdicts (size controlled, power monotone, three-way over-rejecting on a
trend) should agree. Tighter agreement with the reference numbers comes from
raising `MC_N`.

## Notes

- The IPT is **deliberately conservative** under heavy tails and at small
  cluster counts: its p-value lives on the grid `{1, …, K+1}/(K+1)`, so the
  attainable size at `alpha = 0.05` can sit well below 5% when `K` is small
  (e.g. `K = n − 1` with `n = 30` makes `1/(K+1) ≈ 0.033` the only value that
  rejects). Under-rejection is the safe direction; this is a feature, not a
  miscalibration.
- Results are deterministic given the seeds (base seed 20260713; sim `s` draws
  data at `20260713 + s` and fits at `100000 + s`). Exact bit-reproducibility
  additionally depends on `RNGkind()` (the R default is assumed).
- The DGPs are the paper's: the random-feature error model and the Section 7.1
  dyadic design of Guo, Toulis & Wang (2026), plus panel (arbitrary common
  time trend, condition InvB) and fully exchangeable three-way (InvA)
  variants. See `mc_lib.R`.
