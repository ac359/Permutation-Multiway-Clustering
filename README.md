# mwperm

**Finite-sample-exact tests and confidence intervals for regression under
multi-way clustering, panels, replicated layouts, and missing cells.**

`mwperm` implements the invariant permutation test (IPT) of Guo, Toulis & Wang
(2026), built on the residual permutation test of Wen, Wang & Wang (2025). It
is pure base R, with no compiled code.

## Why use this instead of clustered standard errors?

With few clusters, multi-way cluster-robust standard errors (`sandwich::vcovCL`)
and the wild cluster bootstrap rely on asymptotics *in the number of clusters*.
When you have 20–40 countries, firms, or schools per dimension, that
approximation can fail badly, and it fails in the dangerous direction: the test
over-rejects, so you find effects that are not there.

`mwperm` does not approximate. Its p-value is **exact in finite samples** —
the guarantee holds at any number of clusters, and comes from the proof of
Theorem 1 in Guo et al. (2026) rather than from a limit.

The price is resolution, not validity: with `K + 1` permutations the p-value
can only take values `1/(K+1), 2/(K+1), …, 1`. See
[Resolution](#resolution-how-many-clusters-do-you-need).

In this package's own Monte Carlo at α = 0.05, with a single permutation group
per fit (`n_reps = 1`, the configuration Theorem 1 covers), empirical size came
in at or slightly below nominal — the safe direction — across every design
tested: 0.044 dyadic (25 × 25, K = 24, 1000 replications), 0.046 panel
(25 × 25 × 6, K = 24, 500), 0.049 three-way (21 × 21 × 21, K = 20, 900), and
0.035 under heavy-tailed errors (25 × 25 dyadic, K = 24, 600); Monte-Carlo
standard errors 0.006–0.009. A naive OLS test on the same dyadic design
rejected **43.8%** of the time at a nominal 5%. The shipped script
`inst/replication/06_size_by_design.R` covers the remaining designs, at 1000
replications and `n_reps = 1` each: 0.048 replicated layout (10 × 10 cells,
K = 21), 0.043 irregular layout (with a cell-constant covariate, K = 24, and
with one that varies within cells, K = 21), 0.048 incomplete array (40 × 40
minus its diagonal, K = 19) and 0.009 incomplete panel (26 × 26 × 4, K = 22).
Other simulations quoted for this package (the rest of `inst/replication/`,
the software paper) use different sizes, error laws and `n_reps`, and report
different numbers for the same designs; each states its own `n`, `K` and
`n_reps`.

## What it assumes

One assumption per test, and it is about the **errors**, not the covariates.
For every permutation test the error array must be **exchangeable** under
relabelling of the clustering dimensions, conditional on the covariates (Guo
et al., 2026, Assumption 1).

Informally: if you shuffle the country labels, the joint distribution of the
unobserved errors should not change.

This holds, for example, under a two-way random-effects structure

$$\varepsilon_{ij} = \eta_i + \xi_j + u_{ij}$$

with $\eta_i,\xi_j,u_{ij}$ i.i.d. given the covariates. It does **not** require
normality or independence across cells. It does require that the error
variance not depend on the cluster identity or on the covariates: relabelling
the clusters relabels the variance pattern, so covariate-dependent
heteroskedasticity breaks exchangeability *given the covariates* and makes
the permutation test over-reject. For that case the package has a second
invariance route, `mwperm_dyadic_het()`, which replaces relabelling by joint
row-and-column **sign changes** — a sign flip changes no variance, so any
heteroskedasticity is fine, but the errors must then be **symmetric about
zero**. Neither assumption implies the other; see
[Extensions](#extensions).

Critically, neither assumption restricts the **covariate distribution**.
Covariates may be sparse, irregular, or heavy-tailed — exactly the cases
where cluster-robust standard errors are least reliable.

The assumption you must think hardest about is whether a dimension is really
exchangeable. **Time usually is not** — errors are autocorrelated over time, so
permuting periods is invalid. That is what `mwperm_panel()` exists for; see
[Choosing the right design](#choosing-the-right-design).

## Installation

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("ac359/Permutation-Multiway-Clustering")
```

or from a local clone of this repository:

```r
remotes::install_local("mwperm")        # or, in a shell: R CMD INSTALL mwperm
```

Requires R >= 3.6.0. Imports only `stats`, `graphics`, `grDevices` and
`parallel` — all base R. `knitr` and `rmarkdown` (vignette) and `RhpcBLASctl`
(BLAS thread pinning under `n_cores`) are optional.

The **test suite ships with the package**: `tests/` is tracked in this
repository and included in the built tarball, so a clone or an
`R CMD build`/`R CMD check` runs it. The tests are plain base-R `stopifnot`
scripts (no `testthat` dependency) and are deliberately kept fast enough to run
under `R CMD check`. To run them against an installed copy, from the
package root:

```bash
R CMD INSTALL .
for f in tests/test-*.R; do Rscript "$f"; done
```

The top level of `tests/` holds one file per user-facing function --
`test-dyadic.R` tests `mwperm_dyadic()`, and so on -- with the machinery
underneath in `tests/lower-level-tests/`, which `test-lower-level.R` runs.
[`tests/README.md`](tests/README.md) maps the whole suite.

## Quick start

Hand over the data once; `mwperm()` detects the design and runs the matching
test.

```r
library(mwperm)
data(trade_dyadic)        # synthetic 40 x 40 gravity cross-section

fit <- mwperm(y = "log_trade", d = "log_dist",
              x = c("log_gdp_i", "log_gdp_j"),
              index = c("importer", "exporter"),
              data = trade_dyadic, n_reps = 15, seed = 1)
```
```
Detected design: dyadic (2 indices, one observation per cell, complete array)
  -> running mwperm_dyadic(y, d, x, row = importer, col = exporter)
```
```r
fit
```
```
Invariant permutation test (mwperm)
------------------------------------
Design       : dyadic
Auto-detected: dyadic (2 indices, one observation per cell, complete array)
Clusters     : row=40, col=40 (1600 observations)
Permutations : K = 39  (group order 40, 15 reps)
Resolution   : p-values are multiples of 1/40 = 0.025 per rep; reported floor 0.025

  log_dist     OLS estimate = -0.8985   95% IPT CI [-1.246, -0.5485]

H0: beta = 0    p-value = 0.025
Decision     : reject at alpha = 0.05
```

The true coefficient in this synthetic data is −1.0, which the interval covers.

```r
summary(fit)     # tidy data frame, one row per coefficient
confint(fit)     # inverted-test interval, with a "method" attribute
coef(fit)        # OLS point estimate
nobs(fit)        # observations used
plot(fit)        # OLS estimate against the inverted IPT interval
```

```
      term ols_estimate ols_se_naive ipt_ci_low ipt_ci_high p_value
1 log_dist   -0.8985026   0.08728893  -1.246471  -0.5484962   0.025
```

To see what *would* run — design, index roles, balance, attainable resolution —
without any computation:

```r
mwperm_check(index = c("importer", "exporter"), data = trade_dyadic)
```
```
mwperm design diagnosis
------------------------------
Detected design : dyadic (2 indices, one observation per cell, complete array)
Roles           : row = importer, col = exporter
Dimensions      : importer (40) x exporter (40) | 1600 observations
Balance         : complete
Resolution      : default K = 39, so p-values are multiples of 1/40 = 0.025
                  -> fine enough for a 95% confidence set at alpha = 0.05
Would run       : mwperm_dyadic(y, d, x, row = importer, col = exporter)
  ? design = "dyadic_het" runs the sign-flip test instead, which trades
    exchangeability of the errors for symmetry about zero and so
    tolerates arbitrary heteroskedasticity (see ?mwperm_dyadic_het)
```

There is also a formula interface, and the design-specific functions are fully
supported for direct use:

```r
mwperm_formula(log_trade ~ log_dist | log_gdp_i + log_gdp_j,
               data = trade_dyadic, index = c("importer", "exporter"),
               n_reps = 15, seed = 1)

with(trade_dyadic,
     mwperm_dyadic(y = log_trade, d = log_dist,
                   x   = cbind(log_gdp_i, log_gdp_j),
                   row = importer, col = exporter,
                   n_reps = 15, seed = 1))
```

Dispatching through `mwperm()` and calling the design function directly with
the same seed return **identical** objects.

## Choosing the right design

This is the decision that matters. For a finite-sample-exact method the main
practical risk is not arithmetic — it is testing under an exchangeability
assumption your data does not satisfy.

Every design tests the same null in the same model. What differs is **which
labels are permuted**, and therefore **which errors have to be exchangeable**.
Find your data in the first column, then satisfy yourself that the error model
in the third is one you are willing to assume.

| Your data | Use | Exact whenever the errors can be written | What is permuted |
|---|---|---|---|
| Two crossed dimensions, one observation per cell | `mwperm_dyadic()` | $\varepsilon_{ij} = \eta_i + \xi_j + u_{ij}$ | row and column labels, jointly |
| Three crossed dimensions, **all** exchangeable | `mwperm_threeway()` | $\varepsilon_{ijl} = \eta_i + \xi_j + \zeta_l + u_{ijl}$ | all three index sets, jointly |
| Two dimensions observed **over time** | `mwperm_panel()` | $\varepsilon_{ijt} = \eta_i + \xi_j + \zeta_t + u_{ijt}$, with $\zeta_t$ **arbitrary** | rows and columns — the *same* relabelling in every period; $t$ is never moved |
| Repeated **independent** observations per cell | `mwperm_layout()` | $\varepsilon_{ijl} = \eta_{ij} + u_{ijl}$, with $\eta_{ij}$ **arbitrary** | the replicates $l$ inside each cell, drawn independently per cell |
| Repeated observations per cell, but they are **time periods**, or `d` is constant within a cell | `mwperm_irregular()` | $\varepsilon_{ijl} = \eta_i + \xi_j + \zeta_l + u_{ijl}$, with $\zeta_l$ **arbitrary**, and which cells clear $L_0$ not depending on $y$ | whole cells, across rows and columns; the within-cell slot $l$ is never moved |
| Two dimensions over time, but the array has **holes** | `mwperm_panel_missing()` | $\varepsilon_{ijt} = \eta_i + \xi_j + \zeta_t + u_{ijt}$, with $\zeta_t$ **arbitrary**, and the mask $M \perp\mkern-10mu\perp \varepsilon \mid \mathbf{X}, \mathbf{D}$ | rows and columns within each fully observed block, the same relabelling in every period; $t$ is never moved |
| Two dimensions with **missing cells** | `mwperm_missing()` | $\varepsilon_{ij} = \eta_i + \xi_j + u_{ij}$, and the mask $M \perp\mkern-10mu\perp \varepsilon \mid \mathbf{X}, \mathbf{D}$ | rows and columns *within* each fully observed block |
| Two crossed dimensions, one observation per cell, **heteroskedastic** errors | `mwperm_dyadic_het()` | $\varepsilon_{ij} = \sigma_{ij}\, u_{ij}$ with $\sigma_{ij}$ **arbitrary** (it may depend on $i$, $j$ and the covariates) and $u_{ij}$ independent and **symmetric** about zero | the signs of row groups and column groups, jointly; nothing is relabelled |
| Not sure | `mwperm()` or `mwperm_check()` | — | detects the structure and tells you |

In every row but the last the $u$ terms are i.i.d. given the covariates — for
`mwperm_layout()`, i.i.d. *within* each cell, which may differ freely from one
another — and each named family of random effects is i.i.d. within itself;
nothing requires normality, and cluster-level variance heterogeneity that is
itself exchangeable (i.i.d. random scales) is fine, but a variance tied to
the cluster identity or the covariates is not — that is the last row's job.
A term marked **arbitrary** is completely unrestricted — no distribution, no
independence, not even randomness — and that freedom is the reason the row
exists. These are *sufficient* models, given as the easiest way to recognise
your setting; the actual requirement is the invariance below, which each of
them implies.

### The invariance each design needs

All eight designs test the same null in the same regression (Guo et al., 2026,
Eq. 12)

$$y_{ijl} = x_{ijl}^{\top}\gamma + d_{ijl}^{\top}\beta + \varepsilon_{ijl},
\qquad i \in [m],\; j \in [n],\; l \in [\ell_{ij}],$$

against $H_0 : \beta = b$, conditional on $(\mathbf{X}, \mathbf{D})$. With
$\ell_{ij} \equiv 1$ this is the dyadic model of
[The model and the null](#the-model-and-the-null); with $\ell_{ij} \equiv \ell$
constant, $l$ is simply a third cluster dimension. What changes from design to
design is only the group of relabellings the test is invariant to.

**`mwperm_dyadic()` — two-way, condition InvA.** For all permutations $\pi$ on
$[m]$ and $\sigma$ on $[n]$,

$$(\varepsilon_{ij}) \overset{d}{=} (\varepsilon_{\pi(i)\sigma(j)})
\mid \mathbf{X}, \mathbf{D}.$$

**`mwperm_threeway()` — three-way, condition InvA (§6.1).** The same, extended
to a third permutation $\psi$ on $[\ell]$:

$$(\varepsilon_{ijl}) \overset{d}{=}
(\varepsilon_{\pi(i)\sigma(j)\psi(l)}) \mid \mathbf{X}, \mathbf{D}.$$

Use it only if the third dimension really is exchangeable — an industry or a
product category, not a year.

**`mwperm_panel()` — condition InvB (§6.2).** Exchangeability is required
*within* each period only, with the **same** $(\pi, \sigma)$ used in every
period:

$$(\varepsilon_{ijt})_{i \in [m], j \in [n]} \overset{d}{=}
(\varepsilon_{\pi(i)\sigma(j)t})_{i \in [m], j \in [n]}
\mid \mathbf{X}, \mathbf{D}.$$

Nothing at all is assumed *across* $t$: the common trend $\zeta_t$ is arbitrary
and the errors may be autocorrelated over time. This is the first
finite-sample-valid test in that setting. `time_fe = TRUE` (the default) adds
period dummies, which are invariant to the within-period permutation and so
cost no validity, to de-bias the point estimate.

**`mwperm_layout()` — within-cell exchangeability (§6.3).** Independently in
each cell $(i, j)$, for permutations $\pi_{ij}$ on $[\ell_{ij}]$,

$$(\varepsilon_{ijl})_{l \in [\ell_{ij}]} \overset{d}{=}
(\varepsilon_{ij\pi_{ij}(l)})_{l \in [\ell_{ij}]} \mid \mathbf{X}, \mathbf{D}.$$

The cell effects $\eta_{ij}$ are unrestricted, so the clusters themselves need
no exchangeability at all — only the replicates inside a cell. The price is
resolution: $K + 1 \le \min_{ij} \ell_{ij}$ unless `L0 =` balances the array.
Note that $l$ must be an *independent replication*. If the $l$-th observation
means the same thing in every cell — a period, a survey wave — then errors are
usually not exchangeable across $l$ at all, and even a benign shared replicate
effect $\zeta_l$ falls outside this argument, because the independent per-cell
permutations change its alignment across cells. Use `mwperm_irregular()`.

**`mwperm_irregular()` — §6.4.** Condition InvB applied blockwise, with the
within-cell level $l$ (the `rep` argument: a period, a wave) playing time's
role. Choose, from the observation pattern alone, a common set $S$ of $L_0$
levels — the $L_0$ levels jointly observed by the most cells — form the mask
$M_{ij} = \mathbf{1}(\text{cell } (i,j) \text{ observes every } l \in S)$,
find disjoint fully observed blocks under it, keep in every retained cell
exactly the observations at the levels in $S$, then require, within each block
$I_q \times J_q$ and with the same $(\pi, \sigma)$ at every level,

$$(\varepsilon_{ijl})_{i \in I_q, j \in J_q} \overset{d}{=}
(\varepsilon_{\pi(i)\sigma(j)l})_{i \in I_q, j \in J_q}
\mid \mathbf{X}, \mathbf{D}.$$

An arbitrary level effect shared across cells is permitted, which is the whole
point: the permutation never moves an observation to a different $l$. (The
paper's printed step (i) masks on the cell count and drops observations at
random, which does not keep $l$ aligned across cells; see `?mwperm_irregular`.)
It also needs the mask condition below: which cells observe $S$ must not depend
on the outcomes.

**`mwperm_panel_missing()` — §6.2 blockwise, under §5.** An incomplete panel:
the mask keeps the pairs observed in *every* period,
$M_{ij} = \mathbf{1}(\text{pair } (i,j) \text{ observed in all } T)$, the
biclique search cuts it into disjoint fully observed blocks, and within each
block $I_q \times J_q$ condition InvB is required with the same
$(\pi, \sigma)$ in every period:

$$(\varepsilon_{ijt})_{i \in I_q, j \in J_q} \overset{d}{=}
(\varepsilon_{\pi(i)\sigma(j)t})_{i \in I_q, j \in J_q}
\mid \mathbf{X}, \mathbf{D}.$$

Plus Assumption 4 on the mask, below. Nothing is assumed across $t$, so the
trend and any serial correlation are free, exactly as in `mwperm_panel()`. On a
complete array the mask is all ones and the construction *is* `mwperm_panel()`'s:
given the same row and column groups the two build identical permutations,
which the test suite checks directly. (They draw those groups at different
sub-seed offsets, so the same `seed` gives different, equally valid, draws.)

**`mwperm_missing()` — §5.** Cell $(i, j)$ is observed iff $M_{ij} = 1$, and
Assumption 4 requires

$$M \perp\mkern-10mu\perp (\varepsilon_{ij})_{i, j} \mid \mathbf{X}, \mathbf{D}.$$

The mask may depend on $\mathbf{X}$ and $\mathbf{D}$ in any way whatsoever, but
not on $y$. Inference then restricts to disjoint fully observed blocks
$I_q \times J_q$ and applies the two-way invariance inside each; cells outside
the blocks are discarded, and the fit reports how many.

**`mwperm_dyadic_het()` — sign symmetry (§2: any invariance group).** For all
row signs $s \in \lbrace -1, +1 \rbrace^m$ and column signs
$t \in \lbrace -1, +1 \rbrace^n$,

$$(\varepsilon_{ij}) \overset{d}{=} (s_i\, t_j\, \varepsilon_{ij})
\mid \mathbf{X}, \mathbf{D}.$$

No label is moved, so no variance is moved: $\sigma_{ij}$ may depend on $i$,
on $j$ and on the covariates in any way. What is assumed instead is symmetry
of each error about zero (conditionally on the covariates), which InvA never
asked for. Neither condition implies the other, and the sign-flip group is
smaller and coarser than a relabelling group, so under homoskedastic
exchangeable errors this test has less power than `mwperm_dyadic()`. It is a
different bet, not a safer one.

### Still unsure? Four questions

1. **Is one of your indices time, or otherwise ordered?** Then it must never be
   permuted. One observation per $(i, j, t)$ and a complete array →
   `mwperm_panel()`; the same with holes in it → `mwperm_panel_missing()`;
   several observations per $(i, j)$ that are really periods →
   `mwperm_irregular()`.
2. **Do you have more than one observation per $(i, j)$ cell?**
   `mwperm_layout()` if they are exchangeable replicates;
   `mwperm_irregular()` if they are periods, or if `d` is constant within a
   cell — within-cell permutation then has *no power*, because residualizing
   removes all of `d`'s variation.
3. **Are cells missing?** `mwperm_missing()` without a time dimension,
   `mwperm_panel_missing()` with one. Neither is an error path; both are
   different, valid procedures that trade discarded cells for exactness.
4. **Otherwise:** `mwperm_dyadic()` for two indices, `mwperm_threeway()` for
   three genuinely exchangeable ones. If the error *variance* plausibly
   depends on the clusters or the covariates and you are willing to assume
   symmetric errors, `mwperm_dyadic_het()` instead — nothing in the data can
   make that call for you, so `mwperm()` never does.

`mwperm_check(index = ...)` answers all four from the data, prints the
diagnosis and the attainable resolution, and computes nothing.

> **The most costly mistake is treating time as an exchangeable third
> dimension.** Running the three-way test on panel data is *invalid* under
> serial correlation. On a trending null design, the three-way test rejected
> 88% of the time at a nominal 20%; the panel test on the same data rejected
> 19.7%. Running the panel test on genuinely three-way data is merely less
> powerful. So when in doubt, prefer `panel`.

## Function map

| Purpose | Function |
|---|---|
| **Auto-detect and dispatch** | `mwperm()` |
| **Diagnose only (no computation)** | `mwperm_check()` |
| Formula interface | `mwperm_formula()` |
| Two-way / dyadic clustering | `mwperm_dyadic()` |
| Three-way clustering | `mwperm_threeway()` |
| Panel (two-way + arbitrary time trend) | `mwperm_panel()` |
| Incomplete or unbalanced panel | `mwperm_panel_missing()` |
| Replicated two-way layout (`L0=` to balance) | `mwperm_layout()` |
| Irregular layout: repeats are periods, or `d` is cell-level | `mwperm_irregular()` |
| Incomplete array (missing cells) | `mwperm_missing()` |
| Two-way / dyadic clustering with **heteroskedastic** errors (sign-flip test, IPT-Het) | `mwperm_dyadic_het()` |
| Permutation-group construction (Algorithm 1) | `build_perm_set()` |
| Sign-flip-group construction (order $2^{n_{\mathrm{flip}}-1}$) | `build_flip_set()` |
| Fully observed biclique finder (greedy/exact) | `find_bicliques()` |
| Journal-dimension figure export | `mwperm_save()` |

Every test function returns an object of class `"mwperm"` with `print()`,
`summary()`, `confint()`, `coef()`, `nobs()` and `plot()` methods.

## How `mwperm()` chooses, and when it warns

Forks the data itself can settle are resolved **silently**:

- two indices, one observation per cell → `dyadic` if the array is complete,
  `missing` if it has holes;
- three indices, or two indices plus `time =` → `panel` if the array is
  complete, `panel_missing` if it has holes;
- two indices with repeated cells, or `rep =` → `layout`.

`mwperm_irregular()` is never chosen automatically: whether the repeats in a
cell are exchangeable replicates (`layout`) or aligned levels such as periods
(`irregular`) is an assumption, not a fact about the data. Pass
`design = "irregular", L0 = ...` — and `rep =` naming the level — to run it.
`mwperm()` does warn when `d` is constant within every cell, since the layout
test then has no power and the irregular one is the design for that case.

`mwperm_dyadic_het()` is never chosen automatically either, for a stronger
reason: heteroskedasticity leaves **no trace in the clustering structure**, so
there is nothing to detect. On a complete dyadic array `mwperm_check()` prints
one line offering `design = "dyadic_het"`; pass it (with `n_flip =` if you
want other than the default 8) to run the sign-flip test.

Two forks depend on an exchangeability assumption the data *cannot* reveal.
There, `mwperm()` defaults to whichever choice stays valid under the wider set
of error processes, **warns**, and tells you how to override:

- **Three complete crossed indices** → treated as a **panel**, with the
  time-like index held fixed. Force `design = "threeway"` only if all three
  dimensions really are exchangeable.
- **Repeated `(i, j)` cells** → treated as within-cell replication
  (**layout**). If the repeats are time periods, pass `time =` instead.

The time role is assigned from an explicit `time =` tag, then a time-like
column *name*, then time-like *values*. A name-only match that the values do
not corroborate raises a warning, because permuting a true time dimension
over-rejects badly.

## Reading the output

The printed block deliberately labels quantities by **provenance**, because two
different methods are involved:

- `OLS estimate` — the ordinary least-squares point estimate. The IPT is a test,
  not an estimator; this is the usual OLS number.
- `IPT CI` — obtained by **inverting the permutation test**, not as
  `estimate ± 1.96 × SE`. It need not be symmetric about the estimate.
- `ols_se_naive` in `summary()` — the naive homoskedastic OLS standard error.
  It is used internally only to centre and scale the interval search. **It is
  not an inferential quantity — do not report it.**
- `Resolution` — the p-value grid (below).

## Resolution: how many clusters do you need?

The p-value is exact but **discrete**. With group order `K + 1` it can only
take multiples of `1/(K+1)`, so the smallest value attainable is `1/(K+1)`.

Two consequences:

1. To reject at level α you need `1/(K+1) ≤ α`. For α = 0.05 that means
   **K + 1 ≥ 20**, i.e. at least **20 levels in the smallest permuted
   dimension**. With fewer, the p-value is still valid — it simply cannot get
   small enough to reject, and `mwperm` says so in a note and returns no
   confidence interval.
2. A p-value exactly equal to `1/(K+1)` is the *strongest evidence the design
   can produce*, not a precise number.

`K` defaults to `min(permuted dimensions) − 1`, capped at 199. `mwperm_check()`
reports the attainable resolution before you run anything.

`aggregate = "median2"` reports `min(1, 2 × median)`, so its smallest
attainable p-value is `2/(K+1)`: it needs `K + 1 ≥ 2/α`, i.e. **40 levels at
α = 0.05**, twice what the default rule needs. Below that it cannot reject, and
the fit says so in a note and returns no confidence set.

`mwperm_dyadic_het()` is the one exception to "K + 1": its group order is
$2^{n_{\mathrm{flip}}-1}$, so the floor is $1/2^{n_{\mathrm{flip}}-1}$ and a 95%
interval needs **`n_flip ≥ 6`** (order 32); the default `n_flip = 8` gives
order 128 and a floor of 0.0078. The notes and `mwperm_check()` say so in
those terms. See [Extensions](#extensions).

## The model and the null

For the dyadic regression model (Guo et al., 2026, Eq. 1)

$$y_{ij} = x_{ij}^{\top}\gamma + d_{ij}^{\top}\beta + \varepsilon_{ij},
\qquad i, j = 1, \dots, n,$$

with outcome $y_{ij}$, covariates of interest $d_{ij}\in\mathbb{R}^{d}$,
nuisance covariates $x_{ij}\in\mathbb{R}^{p}$ and unobserved errors
$\varepsilon_{ij}$, the package tests

$$H_0 : \beta = b$$

exactly, for any fixed $b$ (default $0$), conditional on $(\mathbf{X},
\mathbf{D})$. Assumption 1 is

$$(\varepsilon_{ij})_{i,j\in[n]} \overset{d}{=}
(\varepsilon_{\pi(i)\sigma(j)})_{i,j\in[n]} \mid \mathbf{X}, \mathbf{D}$$

for all permutations $\pi,\sigma$ on $[n]$.

## How it works

Stacked (Guo et al., 2026, Eq. 7) as
$\mathbf{y} = \mathbf{X}\gamma + \mathbf{D}\beta + \boldsymbol{\varepsilon}$
with $N = n^2$ rows, the test has three steps.

**1. Permutation group.** `build_perm_set()` (Algorithm 1) draws a random,
algebraically closed **block-cyclic group** of $K+1$ two-way permutations

$$\mathcal{G} = \lbrace(\pi_0,\sigma_0),\dots,(\pi_K,\sigma_K)\rbrace,
\qquad \pi_0 = \sigma_0 = \mathrm{Id}.$$

Closure under composition is what the exactness proof requires — a set of
independently drawn permutations would not do.

**2. Partialling out.** For each $k$, form $V_k$ spanning the orthogonal
complement of *both* the nuisance design and its permuted copy,

$$V_k^{\top}\mathbf{X} = 0, \qquad V_k^{\top}\mathbf{X}_{\pi_k,\sigma_k} = 0,$$

a Frisch–Waugh–Lovell projection that removes $\gamma$ without assuming
anything about it. Its dimension is
$N - \mathrm{rank}([\mathbf{X} \mid \mathbf{X}_{\pi_k,\sigma_k}])$.
The paper states this as $V_k \in \mathbb{R}^{N \times (N-2p)}$, which is
the full-rank case; here the stack is rank-deficient *by construction*, since
a permutation maps the intercept to itself (and, in a panel with time effects,
maps the period dummies among themselves), so the projection keeps strictly
more than $N-2p$ dimensions. `mwperm` uses the rank, which is what the two
orthogonality conditions actually ask for and is the only well-defined reading
when the stack is rank-deficient. Then compute

$$a_k = \lVert \mathbf{D}^{\top} V_k V_k^{\top} \mathbf{y}\rVert,
\qquad
b_k = \lVert \mathbf{D}^{\top} V_k V_k^{\top}
\mathbf{y}_{\pi_k,\sigma_k}\rVert.$$

**3. Randomization p-value.** (Guo et al., 2026, Eq. 10)

$$\mathrm{pval} = \frac{1}{K+1}\left(1 + \sum_{k=1}^{K}
\mathbb{1}\left\lbrace \min_{1\le j\le K} a_j \le b_k \right\rbrace\right).$$

Because $\mathcal{G}$ is closed, the identity statistic is exchangeable with the
permuted ones under $H_0$, which makes this exact. The minimum over $j$
(*minorization*) is what keeps it valid under heavy tails, at the cost of some
conservatism.

The permutation group is random, so repeated draws are aggregated by the
**median** p-value across `n_reps` repetitions (Guo et al., 2026, Remark 1).
`n_reps` defaults to 10.

## Confidence sets by test inversion

`confint()` returns

$$\lbrace b \in \mathbb{R}^{d} : \mathrm{pval}(b) > \alpha \rbrace,$$

the set of nulls the test does not reject. Coverage is inherited directly from
the validity of the test — no separate argument is needed, and no normal
approximation is used.

For a single coefficient this set is computed **exactly**, not by search. The
p-value is a step function of $b$ whose jumps solve
$|v_k - W_k b| = |u_j - M_j b|$ for known constants read off the cached
projections, so `mwperm` evaluates it at every jump and at one point inside
every interval between jumps. There is no bracketing assumption and no
bisection tolerance. The set need not be connected — its components are
returned in `fit$conf_set`, a two-column matrix of end points — and `confint()`
reports their hull, which is conservative when there is more than one
component. `fit$ci_method` records which route produced the set (`"exact"`,
`"grid"`, or the `"bisection"` fallback used when the exact candidate count
would be prohibitive).

**What the end points mean.** `conf_set` and `conf_int` are the *closure* of
$\lbrace b : \mathrm{pval}(b) > \alpha \rbrace$. Because the p-value is a
step function, an acceptance region usually begins and ends strictly between
two of its jumps, and there is no attained value at the boundary to report; the
exact route reports the bounding jump. **A reported end point may therefore be
a value the test rejects, while every point strictly inside the interval is
accepted.** The convention is conservative — it never omits an accepted value —
and it means a `b` sitting exactly on an end point should not be read as "just
inside". (The `"grid"` and `"bisection"` routes report attained accepted
points instead, to the grid spacing or bisection tolerance.)

With `n_reps > 1` there is one p-value per repetition, and one rule is used
everywhere: the reported p-value is $\mathrm{median}_r \mathrm{pval}_r(b)$ and
the confidence set is
$\lbrace b : \mathrm{median}_r\, \mathrm{pval}_r(b) > \alpha \rbrace$. The
test and the interval therefore cannot disagree in the direction that matters:
no value the test accepts falls outside the reported set (the end points
themselves are the closure, above). Setting
`aggregate = "median2"` replaces the median with $\min(1, 2\times$ median$)$ in
both places, which restores the level-$\alpha$ guarantee for `n_reps > 1` at
the cost of a wider set.

Two further practical points. The **level is fixed at fit time** (`alpha`), so
`confint(fit, level = 0.90)` on a 95% fit is an error rather than a silent
re-derivation. And with several coefficients the result is a **joint** region;
`confint()` then reports its *marginal extent*, which is not the same as
separate per-coefficient intervals.

Measured coverage of nominal 95% intervals, at 600 simulations per cell (`inst/replication/03_ci_coverage.R`): 0.993 dyadic, 0.992 panel. Coverage above nominal is the valid direction — the p-value lives on the discrete grid `{1, …, K+1}/(K+1)`, so the inverted set is conservative by construction.

## Extensions

All extensions reuse the same machinery; only the invariance condition and the
construction of $\mathcal{G}$ change (Guo et al., 2026, §6). The invariance
each one needs is in [Choosing the right design](#choosing-the-right-design);
what follows is what each function actually does with it.

**Three-way clustering** (`mwperm_threeway()`) applies Algorithm 1 three times,
once per index set, and composes the three into a joint group under InvA.

**Panel data** (`mwperm_panel()`) applies Algorithm 1 twice, to $[m]$ and
$[n]$, then uses that *same* row/column relabelling in every period with time
held fixed — so the dyadic test runs within each period and the unknown trend
is never disturbed. This is the first finite-sample-valid test of $\beta = 0$
under (InvB). Time fixed effects (`time_fe = TRUE`, the default) de-bias the
point estimate; they are invariant to the within-period permutation, so they do
not disturb validity.

```r
data(trade_panel)
mwperm(y = "log_trade", d = "fta", x = c("log_gdp_i", "log_gdp_j"),
       index = c("importer", "exporter", "year"),
       data = trade_panel, seed = 1)
```
```
Detected design: panel ('year' identified as time by name)
  -> running mwperm_panel(y, d, x, row = importer, col = exporter, time = year, time_fe = TRUE)
...
Permutations : K = 21  (group order 22, 10 reps)
Resolution   : p-values are multiples of 1/22 = 0.045 per rep; reported floor 0.045
               the median of 10 reps can fall between grid points

  fta          OLS estimate = 0.6774   95% IPT CI [0.442, 0.8803]

H0: beta = 0    p-value = 0.045
Decision     : reject at alpha = 0.05
```

**Replicated two-way layouts** (`mwperm_layout()`) apply Algorithm 1 once per
cell, on $[\ell_{ij}]$, and permute only *within* cells — appropriate when $l$
indexes independent replications. For unbalanced layouts, `L0` keeps cells with
$\ell_{ij}\ge L_0$ and uniformly downsamples each to exactly $L_0$ replicates
(reproducibly, via `seed`). Note that the `L0` threshold itself comes from
Section 6.4 of the paper, not Section 6.3; `mwperm_layout()` uses it to balance
the array and then runs the Section 6.3 within-cell test.

**Incomplete panels** (`mwperm_panel_missing()`) combine the two: the mask
keeps the (i, j) pairs observed in every period, the biclique search cuts it
into disjoint fully observed blocks, and Procedure 2 runs inside each with the
period held fixed. `mwperm_panel()` refuses an incomplete array outright, and
its error says so. What it costs is cells: a pair observed in five of six
years is dropped whole, and a block has to be complete in both margins, so
thinning a handful of pairs can cost a large share of the array — the fit
reports exactly how much. Dropping the sparsest *periods* before calling is
often the better trade.

**Irregular layouts** (`mwperm_irregular()`) are the Section 6.4 procedure,
and cover the two cases where within-cell permutation fails: the replication
index is really *time* (so within-cell permutation is **invalid**), or
$d_{ijl}$ is constant within each cell (so it has **no power**). It chooses a
common set $S$ of $L_0$ within-cell levels from the observation pattern, forms
the mask $M_{ij} = 1\lbrace \text{cell } (i,j) \text{ observes every level in }
S \rbrace$, runs the biclique search on $M$, keeps in each retained cell
exactly the observations at the levels in $S$, and then applies Procedure 2
*across* cells with the level held fixed — cell $(i,j)$ level $l$ maps to cell
$(\pi(i),\sigma(j))$ level $l$, the same device the panel test uses for time.
It therefore needs exchangeability across $(i,j)$ within each level, not
within-cell exchangeability. With `rep = NULL` the levels are the order of
appearance within the cell and the mask is the paper's
$1\lbrace \ell_{ij} \ge L_0 \rbrace$.

**Heteroskedasticity: the sign-flip test** (`mwperm_dyadic_het()`, IPT-Het).
Section 2 of the paper notes that the partialling-out and minorization
argument is not specific to permutations and holds under *any* invariance
group. This test keeps everything — the model, the null, the statistic
$a_k(b) = \lVert D^\top V_k V_k^\top (y - Db) \rVert$ against
$b_k(b) = \lVert D^\top V_k V_k^\top S_k (y - Db) \rVert$, the minorized
p-value, the exact confidence set by inversion — and replaces the permutation
group $\mathcal{G}$ by a group of **sign flips**: each row cluster and each
column cluster is assigned at random to one of $n_{\mathrm{flip}}$ flip groups, a
sign vector $s \in \lbrace -1,+1 \rbrace^{n_{\mathrm{flip}}}$ multiplies cell
$(i,j)$ by $s_{g_1(i)}\, s_{g_2(j)}$, and element $k$ residualises on
$[X \mid S_k X]$ instead of $[X \mid X_{\pi_k\sigma_k}]$.

- **Assumption.** $(\varepsilon_{ij}) \overset{d}{=} (s_i t_j \varepsilon_{ij})
  \mid \mathbf{X}, \mathbf{D}$: symmetry about zero under joint row-and-column
  sign changes. A sign flip changes no variance, so **arbitrary
  heteroskedasticity** — in $i$, in $j$, in the covariates — is permitted;
  skewed errors are not. Exchangeability is *not* required, and this
  assumption does not imply it, nor the reverse.
- **Group order and resolution.** Because the sign enters as a *product*,
  $s$ and $-s$ induce the same transformation: the $2^{n_{\mathrm{flip}}}$ sign
  vectors give only $2^{n_{\mathrm{flip}}-1}$ distinct elements, and
  `build_flip_set()` enumerates each exactly once. The p-value therefore lives
  on multiples of $1/2^{n_{\mathrm{flip}}-1}$, and a 95% confidence set needs
  $2^{n_{\mathrm{flip}}-1} \ge 20$, i.e. **`n_flip ≥ 6`**. The cost is one
  projection per non-identity element, $2^{n_{\mathrm{flip}}-1} - 1$ per
  repetition — exponential in `n_flip` where the permutation designs are
  linear in `K` — so the default is a fixed `n_flip = 8` (order 128, floor
  0.0078), capped at the smaller cluster count, and values above 20 are
  refused.
- **The trade-off.** This is **not a strict upgrade** over `mwperm_dyadic()`.
  On the authors' heteroskedastic gravity design (25 clusters per side, error
  sd increasing in the gravity mean and in distance) the permutation test's
  size rose with the heteroskedasticity, to 0.06–0.09 at nominal 0.05 in this
  package's 300-replication check (0.13 in the authors' own run), while the
  sign-flip test stayed at 0.02 or below. Under homoskedastic errors, at
  $\beta = 0.15$, the permutation test's power was 0.96 against the sign-flip
  test's 0.86 (0.98 against 0.90 in the authors' run): a smaller, coarser
  group buys robustness with power. Use it when you have reason to doubt
  exchangeability and can defend symmetry; otherwise stay with
  `mwperm_dyadic()`.

```r
with(trade_dyadic,
     mwperm_dyadic_het(y = log_trade, d = log_dist,
                       x = cbind(log_gdp_i, log_gdp_j),
                       row = importer, col = exporter,
                       n_flip = 6, n_reps = 3, seed = 1))
```
```
Design       : dyadic (sign-flip / heteroskedasticity-robust)
Clusters     : row=40, col=40 (1600 observations)
Sign flips   : n_flip = 6 flip groups  (group order 2^5 = 32, 3 reps)
Resolution   : p-values are multiples of 1/32 = 0.031 per rep; reported floor 0.031

  log_dist     OLS estimate = -0.8985   95% IPT CI [-1.16, -0.6251]

H0: beta = 0    p-value = 0.031
Decision     : reject at alpha = 0.05
```

`mwperm(..., design = "dyadic_het", n_flip = 6)` and the formula interface
reach the same function and return the identical object.

## Missing cells

With an observation mask $M$ satisfying
$M \perp\mkern-10mu\perp \boldsymbol{\varepsilon} \mid \mathbf{X},\mathbf{D}$
(Assumption 4), `mwperm_missing()` restricts the permutation to disjoint,
fully observed blocks

$$F_M = \lbrace I_q \times J_q \rbrace_{q=1}^{Q},
\qquad I_q \cap I_{q'} = J_q \cap J_{q'} = \varnothing \quad (q \ne q'),$$

which are **bicliques** in the row–column bipartite graph induced by $M$, and
pools the residual statistics across blocks (Procedure 2).

Two things to expect. First, **cells outside the selected blocks are
discarded** — routinely well over half. That is the mechanism, not a
malfunction: it is what leaves the permutation acting on data whose
exchangeability structure is intact. The fit reports exactly how many cells
were kept.

Second, **the smallest selected block caps `K`**, and therefore the resolution.
A fit can be perfectly valid yet unable to reject at α = 0.05 because one small
block set `K = 4`. `mwperm` says so explicitly in a note; raise `min_block` to
stop small blocks from setting `K`.

Maximum-biclique search is NP-hard, so the default is a greedy heuristic
(`block_method = "exact"` runs branch-and-bound with a node budget, falling
back to greedy). **A sub-maximal block costs power, never validity.** Power
exhibits phase transitions in the missingness rate (Guo et al., 2026, §5).

## Reproducibility

Pass `seed =` and results are exactly reproducible. Repetition `r` uses seed
`seed + r − 1`, and all internal randomness is drawn through a save/restore
wrapper, so **a seeded `mwperm` call never disturbs your global RNG stream**.

With `seed = NULL` the permutations come from the ambient RNG and results vary
between runs. Because the p-value depends on a random group, reporting a seed
alongside a p-value is good practice; `n_reps > 1` with the median aggregation
reduces the dependence on any one draw.

## Performance and parallelism

Cost scales as roughly `n_reps × K × N × (p² + pd + d²)` — linear in the number
of repetitions, the group order, and the sample size. Indicative serial timings
on a laptop, `n_reps = 15` with confidence intervals:

| Design | Clusters per dimension | N | Time |
|---|---|---|---|
| dyadic | 100 | 10,000 | 1.2 s |
| dyadic | 200 | 40,000 | 8.7 s |
| panel (5 periods) | 100 | 50,000 | 12 s |
| panel (5 periods) | 200 | 200,000 | 93 s |

Every test function takes `n_cores =`. Results are **identical** to the serial
run — verified by `identical()` at 2, 4 and 8 cores — because every draw is
derived from an explicit seed rather than from worker scheduling.

Parallelism helps when `n_reps > 1` **and** `seed` is set, which is the axis
that scales (measured 1.8×, 3.1×, 4.6× at 2, 4, 8 cores). With `n_reps = 1`
there is little to gain. Confidence intervals are nearly free: the interval
search reuses cached factorizations and adds only a few percent.

## Limitations

- **Resolution, not validity, is the binding constraint** with few clusters.
  Below 20 levels in the smallest permuted dimension you get a valid p-value
  but no 95% interval.
- **The test is conservative under heavy tails.** The minorized statistic
  under-rejects rather than over-rejects. This is deliberate and is the safe
  direction, but it costs power.
- **Heteroskedasticity is a choice you make, not one the package detects.**
  The permutation tests need exchangeability *given the covariates*;
  `mwperm_dyadic_het()` needs symmetric errors. Nothing in the data settles
  which holds, and the sign-flip test is less powerful when both do.
- **Incomplete panels pay in cells.** `mwperm_panel_missing()` keeps only
  the pairs observed in *every* period, and only those inside fully observed
  blocks, so a few thinly observed pairs can cost a large share of the array;
  the fit reports exactly how much. Dropping the sparsest periods first is
  often the better trade.
- **Exact biclique search does not scale** to large dense-complement masks; it
  falls back to greedy with a warning. Validity is unaffected.
- **`mwperm` tests coefficients; it does not estimate them.** Point estimates
  are OLS. There is no IPT estimator or IPT standard error.
- The p-value depends on a random permutation group, so different seeds give
  slightly different values. Use `n_reps > 1` and report your seed.

## Data

Two synthetic data sets ship with the package and are used in the examples:

- `trade_dyadic` — complete 40 × 40 cross-section of bilateral trade;
- `trade_panel` — balanced 22 × 22 × 6 panel with a free-trade-agreement
  treatment.

Both are generated reproducibly by `data-raw/make_data.R` and carry their
data-generating coefficients in `attr(., "true_coef")`, so examples can be
checked against a known truth. **They contain no real trade statistics.**

## Citation and references

If you use this package, please cite the method paper:

- Guo, W., Toulis, P. & Wang, Y. (2026). *Permutation Inference under Multi-way
  Clustering and Missing Data.* arXiv:2601.08610 [stat.ME].
  <doi:10.48550/arXiv.2601.08610>
- Wen, K., Wang, T. & Wang, Y. (2025). *Residual permutation test for
  regression coefficient testing.* Annals of Statistics 53(2), 724–748.
  <doi:10.1214/24-AOS2479>

## License

MIT. See `LICENSE` for the copyright year and holders, and `LICENSE.md` for the
full text.
