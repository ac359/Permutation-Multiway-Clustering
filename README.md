# mwperm

**Invariant permutation tests for multi-way clustered, panel, and missing-data regression.**

`mwperm` implements the finite-sample-valid permutation test of Guo, Toulis &
Wang (2026), which extends the residual permutation test of Wen, Wang & Wang
(2025, *Annals of Statistics*) to multi-way clustered designs. It tests a
regression coefficient

```
$$y_{ij} = x_{ij}'\gamma + d_{ij}'\beta + \epsilon_{ij}, \quad H_0: \beta = b$$
```

where the errors are *separately exchangeable* across clustering dimensions
(e.g. `eps_ij = eta_i + xi_j + u_ij`). Unlike multi-way cluster-robust
standard errors or the wild cluster bootstrap, the test

- is **exact in finite samples** — valid with as few as ~20 clusters per
  dimension, with no asymptotics in the number of clusters;
- makes **no assumption on the covariate distribution** — covariates may be
  sparse, irregular, or heavy-tailed;
- extends cleanly to **three-way clustering**, **panel data** with an
  arbitrary common time trend, **replicated two-way layouts**, and
  **incomplete arrays** (missing cells).

## Installation

From a local clone of this repository:

```r
# install.packages("remotes")
remotes::install_local("mwperm")        # or: R CMD INSTALL mwperm
```

The package depends only on base R (>= 3.6) and the recommended packages
`stats` and `graphics`.

## Quick start

```r
library(mwperm)
data(trade_dyadic)        # synthetic 40 x 40 gravity panel

fit <- with(trade_dyadic,
            mwperm_dyadic(y = log_trade, d = log_dist,
                          x   = cbind(log_gdp_i, log_gdp_j),
                          row = importer, col = exporter,
                          n_reps = 15, seed = 1))
fit
#>   log_dist   estimate = -0.89   95% CI [-1.25, -0.56]
#>   H0: beta = 0    p-value = 0.025   (reject at alpha = 0.05)

summary(fit)              # tidy data frame
confint(fit)              # inverted-test confidence interval
plot(fit)                 # stability of the randomised p-value across reps
```

## Function map

| Design                                   | Function            |
|------------------------------------------|---------------------|
| Two-way / dyadic clustering              | `mwperm_dyadic()`   |
| Three-way clustering                     | `mwperm_threeway()` |
| Panel (two-way + arbitrary time trend)   | `mwperm_panel()`    |
| Replicated two-way layout (unequal cells; `L0=` to balance) | `mwperm_layout()` |
| Incomplete array (missing cells)         | `mwperm_missing()`  |
| Permutation-group construction (Alg. 1)  | `build_perm_set()`  |
| Fully observed biclique finder (greedy/exact) | `find_bicliques()` |

All test functions return an object of class `"mwperm"` with `print()`,
`summary()`, `confint()` and `plot()` methods.

## How it works

For each member `k` of a randomly drawn, algebraically closed permutation
group (a block-cyclic group of order `K + 1`, built by `build_perm_set()`),
the test forms the residual statistic `||D' V_k V_k' y||`, where the columns
of `V_k` span the orthogonal complement of both the nuisance design `X` and
its permuted copy. Because the group is closed under composition, the
statistic for the identity is exchangeable with the others under `H0`, giving
an exact p-value

```
p = (1 / (K + 1)) * (1 + #{ k : min_j a_j <= b_k }).
```

The smallest attainable p-value is `1 / (K + 1)`, so a 95% confidence
interval (obtained by inverting the test) requires at least 20 clusters per
dimension. Randomised runs are aggregated by the **median** p-value
(`n_reps`).

- **Panel data** (`mwperm_panel()`) applies the *same* row/column permutation
  in every period, holding an arbitrary time trend fixed (condition InvB).
  This is the first finite-sample-valid test for `beta = 0` in such panels.
  Time fixed effects (`time_fe = TRUE`, the default) de-bias the point
  estimate and sharpen the interval.
- **Missing cells** (`mwperm_missing()`) restrict the permutation to disjoint,
  fully observed bicliques and pool the residual statistics across them.

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

- Guo, F. R., Toulis, P. & Wang, Y. (2026). *Permutation inference under
  multi-way clustering and missing data.*
- Wen, K., Wang, T. & Wang, Y. (2025). Residual permutation test for
  regression coefficient testing. *The Annals of Statistics* **53**(2),
  724–748. doi:10.1214/24-AOS2360.

## License

MIT (see `LICENSE`).
