# mwperm

**Invariant permutation tests for multi-way clustered, panel, and missing-data regression.**

`mwperm` implements the finite-sample-valid *invariant permutation test* (IPT) of
Guo, Toulis & Wang (2026), which extends the residual permutation test of Wen,
Wang & Wang (2025, *Annals of Statistics*) to multi-way clustered designs.

For the **dyadic regression model** (Guo et al., 2026, Eq. 1)

$$
y_{ij} = x_{ij}^{\top}\gamma + d_{ij}^{\top}\beta + \varepsilon_{ij},
\qquad i, j = 1, \dots, n,
$$

with outcomes $y_{ij}\in\mathbb{R}$, covariates of interest $d_{ij}\in\mathbb{R}^{d}$,
auxiliary (nuisance) covariates $x_{ij}\in\mathbb{R}^{p}$, and unobserved errors
$\varepsilon_{ij}$, the package tests the significance null

$$
H_0 : \beta = 0
$$

(and, by test inversion, any $H_0^{b}:\beta = b$), where $\beta$ is the parameter of
interest and $\gamma$ is a nuisance parameter.

The test is exact whenever the errors are **conditionally (doubly) exchangeable**
given the covariates — Assumption 1 of Guo et al. (2026):

$$
(\varepsilon_{ij})_{i,j\in[n]}  \overset{d}{=} 
(\varepsilon_{\pi(i)\sigma(j)})_{i,j\in[n]}  \mid  \mathbf{X}, \mathbf{D}
\qquad \text{for all permutations } \pi, \sigma \text{ on } [n] := \lbrace1,\dots,n\rbrace.
$$

This holds, for example, under a two-way random-effects structure
(Guo et al., 2026, Eq. 8)

$$
\varepsilon_{ij} = \eta_i + \xi_j + u_{ij},
$$

with $\eta_i, \xi_j, u_{ij}$ i.i.d. conditional on the covariates. Unlike the
*separate-exchangeability* framework that underlies multi-way cluster-robust
standard errors and the wild cluster bootstrap, Assumption 1 imposes **no restriction
on the covariate distribution** (covariates may be sparse, irregular, or heavy-tailed)
and holds in the *finite sample* rather than the population. Relative to those methods,
the test:

- is **exact in finite samples** — valid with as few as ~20 clusters per dimension,
  with no asymptotics in the number of clusters (Guo et al., 2026, Theorem 1);
- makes **no assumption on the covariate distribution** — covariates may be sparse,
  irregular, or heavy-tailed;
- extends cleanly to **three-way clustering**, **panel data** with an arbitrary common
  time trend, **replicated two-way layouts**, and **incomplete arrays** (missing cells);
- **detects the design for you**: `mwperm()` inspects the clustering structure and
  dispatches to the right test, and `mwperm_check()` prints the diagnosis without
  running anything.

## Installation

From a local clone of this repository:

```r
# install.packages("remotes")
remotes::install_local("mwperm")        # or: R CMD INSTALL mwperm
```

The package depends only on base R (>= 3.6): `stats`, `graphics`, and
`parallel` (for the opt-in `n_cores` argument). `knitr`, `rmarkdown` and
`RhpcBLASctl` are optional (vignette, BLAS pinning).

## Quick start

Hand the data over once; `mwperm()` figures out the design and runs the
matching test:

```r
library(mwperm)
data(trade_dyadic)        # synthetic 40 x 40 gravity cross-section

fit <- mwperm(y = "log_trade", d = "log_dist",
              x = c("log_gdp_i", "log_gdp_j"),
              index = c("importer", "exporter"),
              data = trade_dyadic, n_reps = 15, seed = 1)
#> Detected design: dyadic (2 indices, one observation per cell, complete array)
#>   -> running mwperm_dyadic(y, d, x, row = importer, col = exporter)
fit
#>   log_dist   estimate = -0.89   95% CI [-1.25, -0.56]
#>   H0: beta = 0    p-value = 0.025   (reject at alpha = 0.05)

summary(fit)              # tidy data frame
confint(fit)              # inverted-test confidence interval
plot(fit)                 # stability of the randomised p-value across reps
```

To see what would run — detected design, index roles, balance, attainable
p-value resolution — without any computation:

```r
mwperm_check(index = c("importer", "exporter"), data = trade_dyadic)
```

The design-specific functions remain fully supported for direct use:

```r
fit <- with(trade_dyadic,
            mwperm_dyadic(y = log_trade, d = log_dist,
                          x   = cbind(log_gdp_i, log_gdp_j),
                          row = importer, col = exporter,
                          n_reps = 15, seed = 1))
```

For large problems (thousands of cells, large `K`), every test function takes
`n_cores =` to parallelise the permutation computations; the result is
**identical** to the serial one (every draw is derived from explicit seeds).

## Function map

| Design                                   | Function            |
|------------------------------------------|---------------------|
| **Auto-detect and dispatch**             | `mwperm()`          |
| **Diagnose only (no computation)**       | `mwperm_check()`    |
| Two-way / dyadic clustering              | `mwperm_dyadic()`   |
| Three-way clustering                     | `mwperm_threeway()` |
| Panel (two-way + arbitrary time trend)   | `mwperm_panel()`    |
| Replicated two-way layout (unequal cells; `L0=` to balance) | `mwperm_layout()` |
| Incomplete array (missing cells)         | `mwperm_missing()`  |
| Permutation-group construction (Alg. 1)  | `build_perm_set()`  |
| Fully observed biclique finder (greedy/exact) | `find_bicliques()` |

All test functions return an object of class `"mwperm"` with `print()`,
`summary()`, `confint()` and `plot()` methods.

### How `mwperm()` chooses

Forks the data itself can resolve (complete vs incomplete arrays, repeated
cells) are decided silently. Two forks hinge on an exchangeability assumption
the data *cannot* reveal, and there `mwperm()` defaults to the choice that
stays valid under the widest set of error processes and tells you how to
override: a complete 3-index array is treated as a **panel** (time held
fixed) unless you force `design = "threeway"` — the panel test remains valid
when all three dimensions are exchangeable, while the converse fails — and
repeated `(i, j)` cells are treated as within-cell replication (**layout**)
with a warning, since they could equally be a time dimension you forgot to
pass (`time =` if so).

## How it works

Write the stacked dyadic model (Guo et al., 2026, Eq. 7) as

$$
\mathbf{y} = \mathbf{X}\gamma + \mathbf{D}\beta + \boldsymbol{\varepsilon},
\qquad
\mathbf{y},\boldsymbol{\varepsilon}\in\mathbb{R}^{N},\quad
\mathbf{X}\in\mathbb{R}^{N\times p},\quad
\mathbf{D}\in\mathbb{R}^{N\times d},\quad N = n^2,
$$

where the entries are stacked in lexicographic order, so that cell $(i,j)$ occupies
row $(i-1)n + j$. The test then follows three steps.

**1. Permutation group.** `build_perm_set()` (Algorithm 1) draws a random,
algebraically closed **block-cyclic group** of $K+1$ two-way permutations

$$
\mathcal{G} = \lbrace(\pi_0,\sigma_0), (\pi_1,\sigma_1), \dots, (\pi_K,\sigma_K)\rbrace,
\qquad \pi_0 = \sigma_0 = \mathrm{Id},
$$

closed under composition in each dimension.

**2. Partialling-out and the residual statistic.** For each $k = 1,\dots,K$, the test
forms an orthonormal $V_k\in\mathbb{R}^{N\times(N-2p)}$ whose columns span the
orthogonal complement of *both* the nuisance design and its permuted copy,

$$
V_k^{\top}\mathbf{X} = 0,
\qquad V_k^{\top}\mathbf{X}_{\pi_k,\sigma_k} = 0
$$

(a Frisch–Waugh–Lovell-style projection), and computes the residual statistics

$$
a_k = \lVert \mathbf{D}^{\top} V_k V_k^{\top} \mathbf{y}\rVert,
\qquad
b_k = \lVert \mathbf{D}^{\top} V_k V_k^{\top} \mathbf{y}_{\pi_k,\sigma_k}\rVert,
$$

where $\mathbf{y}_{\pi_k,\sigma_k}$ is $\mathbf{y}$ after permuting rows by $\pi_k$ and
columns by $\sigma_k$. Because $\mathcal{G}$ is closed under composition, the identity
statistic is exchangeable with the permuted ones under $H_0$.

**3. Randomization p-value.** The exact p-value (Guo et al., 2026, Eq. 10) is

$$
\mathrm{pval} = \frac{1}{K+1}
\left(1 + \sum_{k=1}^{K}
\mathbb{1}\left\lbrace \min_{1\le j\le K} a_j \le b_k \right\rbrace\right),
$$

and `confint()` inverts the test to return the $100(1-\alpha)\%$ confidence region

$$
\mathrm{CI} = \lbrace b \in \mathbb{R}^{d} : \mathrm{pval}(b) > \alpha \rbrace.
$$

The smallest attainable p-value is $1/(K+1)$, and Algorithm 1 requires $n \ge K+1$;
hence a 95% confidence interval needs $K + 1 \ge 20$, i.e. **at least 20 clusters per
dimension**. The permutation group is random, so randomized runs are aggregated by the
**median** p-value across `n_reps` repetitions (Guo et al., 2026, Remark 1).

## Extensions

All extensions reuse the same partialling-out / minorization machinery; only the
invariance condition and the construction of $\mathcal{G}$ change (Guo et al., 2026, §6).

**Three-way clustering** (`mwperm_threeway()`). For

$$
y_{ijl} = x_{ijl}^{\top}\gamma + d_{ijl}^{\top}\beta + \varepsilon_{ijl},
\qquad i\in[m],  j\in[n],  l\in[\ell],
$$

under the three-way exchangeability condition (InvA),

$$
(\varepsilon_{ijl})  \overset{d}{=} 
(\varepsilon_{\pi(i)\sigma(j)\psi(l)})  \mid  \mathbf{X}, \mathbf{D},
$$

Algorithm 1 is applied three times. (InvA) holds, e.g., under the full random-effects
model $\varepsilon_{ijl} = \eta_i + \xi_j + \zeta_l + u_{ijl}$.

**Panel data** (`mwperm_panel()`). With a time index $t$,

$$
y_{ijt} = x_{ijt}^{\top}\gamma + d_{ijt}^{\top}\beta + \varepsilon_{ijt},
$$

full three-way exchangeability is implausible because errors are autocorrelated over
$t$. The test instead requires exchangeability across the first two dimensions only
(InvB),

$$
(\varepsilon_{ijt})_{i\in[m],j\in[n]}  \overset{d}{=} 
(\varepsilon_{\pi(i)\sigma(j)t})_{i\in[m],j\in[n]}  \mid  \mathbf{X}, \mathbf{D},
$$

which holds under $\varepsilon_{ijt} = \eta_i + \xi_j + \zeta_t + u_{ijt}$ with $\zeta_t$
an **arbitrary common time trend**. The *same* row/column permutation is applied in every
period, so the procedure runs the baseline dyadic test within each $t$ and holds the time
trend fixed. This is the first finite-sample-valid test of $\beta = 0$ under (InvB). Time
fixed effects (`time_fe = TRUE`, the default) de-bias the point estimate and sharpen the
interval.

**Replicated two-way layouts** (`mwperm_layout()`). When the number of observations per
cell $\ell_{ij}$ varies, permutation is performed only within each cell $(i,j)$ over
$[\ell_{ij}]$. This is valid under the relaxed structure
$\varepsilon_{ijl} = \eta_{ij} + \zeta_l + u_{ijl}$ ($\eta_{ij}$ arbitrary, $\zeta_l$
i.i.d.) — appropriate when $l$ indexes independent replications, as in two-way layouts
from randomized experiments. For unbalanced layouts, the threshold argument `L0` keeps
the cells with $\ell_{ij}\ge L_0$, uniformly downsamples each to exactly $L_0$
replicates (reproducibly via `seed`), and runs the same within-cell test on the
resulting balanced array; tune `L0` (e.g. by grid search) to trade off the number of
eligible cells against within-cell sample size.

**Missing cells** (`mwperm_missing()`). With an observation mask
$M\in\lbrace0,1\rbrace^{n\times n}$ ($M_{ij} = 1$ iff cell $(i,j)$ is observed) satisfying
$M \perp\mkern-10mu\perp \boldsymbol{\varepsilon} \mid \mathbf{X},\mathbf{D}$
(Assumption 4), `find_bicliques()` restricts the permutation to disjoint, fully observed
blocks

$$
F_M = \lbrace I_q \times J_q \rbrace_{q=1}^{Q},
\qquad I_q \cap I_{q'} = J_q \cap J_{q'} = \varnothing \quad (q \ne q'),
$$

which are exactly **bicliques** (complete bipartite subgraphs) in the row–column
bipartite graph induced by $M$. The residual statistics are pooled across blocks. Power
exhibits phase transitions in the missingness rate under Erdős–Rényi-type masking
(Guo et al., 2026, §5).

## Data

Two synthetic data sets ship with the package and are used throughout the
examples and the technical report:

- `trade_dyadic` — a complete 40 x 40 cross-section of bilateral trade;
- `trade_panel` — a balanced 22 x 22 x 6 longitudinal version with a policy
  (free-trade-agreement) treatment.

Both are generated reproducibly by `data-raw/make_data.R` and carry the
data-generating coefficients in `attr(., "true_coef")`. They contain no real
trade statistics.

## References

- Guo, W., Toulis, P. & Wang, Y. (2026). *Permutation Inference under Multi-way
  Clustering and Missing Data.* arXiv:2601.08610 [stat.ME].
- Wen, K., Wang, T. & Wang, Y. (2025). Residual permutation test for regression
  coefficient testing. *The Annals of Statistics* **53**(2), 724–748.
  doi:10.1214/24-AOS2360.
