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

In this package's own Monte Carlo at α = 0.05, empirical size came in at or
slightly below nominal — the safe direction — across every design tested:
0.044 dyadic, 0.046 panel, 0.049 three-way, and 0.035 under heavy-tailed
errors (Monte-Carlo standard errors 0.006–0.009; 500–1000 replications per
cell). A naive OLS test on the same dyadic design rejected **43.8%** of the
time at a nominal 5%.

## What it assumes

One assumption, and it is about the **errors**, not the covariates: the error
array must be **exchangeable** under relabelling of the clustering dimensions,
conditional on the covariates (Guo et al., 2026, Assumption 1).

Informally: if you shuffle the country labels, the joint distribution of the
unobserved errors should not change.

This holds, for example, under a two-way random-effects structure

$$\varepsilon_{ij} = \eta_i + \xi_j + u_{ij}$$

with $\eta_i,\xi_j,u_{ij}$ i.i.d. given the covariates. It does **not** require
homoskedasticity, normality, or independence across cells.

Critically, it imposes **no restriction whatsoever on the covariate
distribution**. Covariates may be sparse, irregular, or heavy-tailed — exactly
the cases where cluster-robust standard errors are least reliable.

The assumption you must think hardest about is whether a dimension is really
exchangeable. **Time usually is not** — errors are autocorrelated over time, so
permuting periods is invalid. That is what `mwperm_panel()` exists for; see
[Choosing the right design](#choosing-the-right-design).

## Installation

From a local clone of this repository:

```r
# install.packages("remotes")
remotes::install_local("mwperm")        # or: R CMD INSTALL mwperm
```

Requires R >= 3.6.0. Imports only `stats`, `graphics`, `grDevices` and
`parallel` — all base R. `knitr` and `rmarkdown` (vignette) and `RhpcBLASctl`
(BLAS thread pinning under `n_cores`) are optional.

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
Resolution   : p-values are multiples of 1/40 = 0.025

  log_dist     OLS estimate = -0.8985   95% IPT CI [-1.246, -0.5486]

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
1 log_dist   -0.8985026   0.08728893  -1.246414  -0.5485613   0.025
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
                  -> fine enough for a 95% confidence set
Would run       : mwperm_dyadic(y, d, x, row = importer, col = exporter)
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

| Your data | Use | Because |
|---|---|---|
| Two crossed dimensions, one observation per cell | `mwperm_dyadic()` | errors exchangeable in both margins |
| Three crossed dimensions, **all** exchangeable | `mwperm_threeway()` | permutes all three |
| Two dimensions observed **over time** | `mwperm_panel()` | holds time fixed; valid under an *arbitrary* common time trend |
| Repeated **independent** observations per cell | `mwperm_layout()` | permutes replicates within cells |
| Two dimensions with **missing cells** | `mwperm_missing()` | restricts to fully observed blocks |
| Not sure | `mwperm()` or `mwperm_check()` | detects and tells you |

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
| Replicated two-way layout (`L0=` to balance) | `mwperm_layout()` |
| Incomplete array (missing cells) | `mwperm_missing()` |
| Permutation-group construction (Algorithm 1) | `build_perm_set()` |
| Fully observed biclique finder (greedy/exact) | `find_bicliques()` |
| Journal-dimension figure export | `mwperm_save()` |

Every test function returns an object of class `"mwperm"` with `print()`,
`summary()`, `confint()`, `coef()`, `nobs()` and `plot()` methods.

## How `mwperm()` chooses, and when it warns

Forks the data itself can settle — complete vs incomplete arrays, repeated
cells, number of indices — are resolved **silently**.

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
anything about it, then compute

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

Two practical points. The **level is fixed at fit time** (`alpha`), so
`confint(fit, level = 0.90)` on a 95% fit is an error rather than a silent
re-derivation. And with several coefficients the result is a **joint** region;
`confint()` then reports its *marginal extent*, which is not the same as
separate per-coefficient intervals.

Measured coverage of nominal 95% intervals: 0.950 dyadic, 0.963 panel.

## Extensions

All extensions reuse the same machinery; only the invariance condition and the
construction of $\mathcal{G}$ change (Guo et al., 2026, §6).

**Three-way clustering** (`mwperm_threeway()`) applies Algorithm 1 three times
under three-way exchangeability (InvA), which holds e.g. under
$\varepsilon_{ijl} = \eta_i + \xi_j + \zeta_l + u_{ijl}$.

**Panel data** (`mwperm_panel()`) requires exchangeability across the first two
dimensions only (InvB):

$$(\varepsilon_{ijt})_{i\in[m],j\in[n]} \overset{d}{=}
(\varepsilon_{\pi(i)\sigma(j)t})_{i\in[m],j\in[n]} \mid \mathbf{X},\mathbf{D},$$

which holds under $\varepsilon_{ijt} = \eta_i + \xi_j + \zeta_t + u_{ijt}$ with
$\zeta_t$ an **arbitrary common time trend**. The *same* row/column permutation
is applied in every period and time is held fixed. This is the first
finite-sample-valid test of $\beta = 0$ under (InvB). Time fixed effects
(`time_fe = TRUE`, the default) de-bias the point estimate; they are invariant
to the within-period permutation, so they do not disturb validity.

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
Resolution   : p-values are multiples of 1/22 = 0.045

  fta          OLS estimate = 0.6774   95% IPT CI [0.4441, 0.8776]

H0: beta = 0    p-value = 0.045
Decision     : reject at alpha = 0.05
```

**Replicated two-way layouts** (`mwperm_layout()`) permute only *within* each
cell $(i,j)$ over $[\ell_{ij}]$, valid under
$\varepsilon_{ijl} = \eta_{ij} + \zeta_l + u_{ijl}$ with $\eta_{ij}$ arbitrary —
appropriate when $l$ indexes independent replications. For unbalanced layouts,
`L0` keeps cells with $\ell_{ij}\ge L_0$ and uniformly downsamples each to
exactly $L_0$ replicates (reproducibly, via `seed`).

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
- **Incomplete or unbalanced panels are not implemented.** Panels must be
  complete balanced arrays. Combining the panel condition with the biclique
  restriction is a natural extension but is not shipped.
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
