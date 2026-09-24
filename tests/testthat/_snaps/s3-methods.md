# print() output is stable (snapshot)

    Code
      print(fit)
    Output
      
      Invariant permutation test (mwperm)
      ------------------------------------
      Design       : dyadic
      Clusters     : row=20, col=20 (400 observations)
      Permutations : K = 19  (group order 20, 5 reps)
      Resolution   : p-values are multiples of 1/20 = 0.050 per rep; reported floor 0.050
      
        fx$d         OLS estimate = 0.3162   95% IPT CI [0.1129, 0.5042]
      
      H0: beta = 0    p-value = 0.050
      Decision     : reject at alpha = 0.05
      

---

    Code
      print(mwperm_missing(md$y, md$d, x = md$x, row = md$i, col = md$j, n_reps = 3,
      seed = 1))
    Output
      
      Invariant permutation test (mwperm)
      ------------------------------------
      Design       : missing (bicliques)
      Clusters     : row=20, col=20 (200 observations)
      Permutations : K = 9  (group order 10, 3 reps)
      Resolution   : p-values are multiples of 1/10 = 0.100 per rep; reported floor 0.100
      
        md$d         OLS estimate = 0.3937
      
      H0: beta = 0    p-value = 0.100
      Decision     : do not reject at alpha = 0.05
      
      Notes:
        - Kept 200 of 380 observed cells (52.6%): the 2 fully observed
          rectangular blocks Procedure 2 could extract. The other 180 cells
          are discarded -- under missingness that loss is what buys exact
          validity, and a larger block would only add power.
        - Block sizes (rows x cols): 10x10, 10x10.
        - Resolution here is set by the smallest selected block: its permuted
          side is 10, so K = 9. Raise `min_block` so that small blocks cannot
          set K -- a higher floor discards more cells but lifts the
          attainable resolution.
        - No 95% confidence interval: the smallest attainable p-value is
          1/(K+1) = 0.1, which is above alpha = 0.05, so no value could be
          excluded and the set would be the whole line. A 95% set needs K + 1
          >= 20 -- that is, at least 20 levels in the smallest permuted
          dimension. The p-value reported above is unaffected and remains
          exact.
      

---

    Code
      print(mwperm_dyadic_het(fx$y, fx$d, x = fx$x, row = fx$i, col = fx$j, n_flip = 6,
      n_reps = 3, seed = 1))
    Output
      
      Invariant sign-flip test (mwperm)
      ------------------------------------
      Design       : dyadic (sign-flip / heteroskedasticity-robust)
      Clusters     : row=20, col=20 (400 observations)
      Sign flips   : n_flip = 6 flip groups  (group order 2^5 = 32, 3 reps)
      Resolution   : p-values are multiples of 1/32 = 0.031 per rep; reported floor 0.031
      
        fx$d         OLS estimate = 0.3162   95% IPT CI [-0.05024, 0.6152]
      
      H0: beta = 0    p-value = 0.062
      Decision     : do not reject at alpha = 0.05
      

