# mwperm_check() prints a diagnosis for every design

    Code
      print(mwperm_check(index = dy[c("i", "j")]))
    Output
      
      mwperm design diagnosis
      ------------------------------
      Detected design : dyadic (2 indices, one observation per cell, complete array)
      Roles           : row = i, col = j
      Dimensions      : i (20) x j (20) | 400 observations
      Balance         : complete
      Resolution      : default K = 19, so p-values are multiples of 1/20 = 0.050
                        -> fine enough for a 95% confidence set at alpha = 0.05
      Would run       : mwperm_dyadic(y, d, x, row = i, col = j)
        ? design = "dyadic_het" runs the sign-flip test instead: valid under
          arbitrary heteroskedasticity for errors that are independent across
          cells and symmetric about zero, but NOT under additive cluster
          effects eta_i + xi_j (see ?mwperm_dyadic_het)
      

---

    Code
      print(mwperm_check(index = md[c("i", "j")]))
    Output
      
      mwperm design diagnosis
      ------------------------------
      Detected design : missing (2 indices, 338 of 400 cells observed)
      Roles           : row = i, col = j
      Dimensions      : i (20) x j (20) | 338 observations
      Balance         : incomplete (338 of 400 cells)
      Resolution      : set by the biclique blocks, so not known until
                        they are found (see find_bicliques)
      Would run       : mwperm_missing(y, d, x, row = i, col = j, min_block = ...)
        - The permutation-group order under missingness is set by the fully
          observed biclique blocks; see find_bicliques() for the achievable
          K.
        ? design = "dyadic_het" runs the sign-flip test on every observed
          cell, with no biclique search and nothing discarded: valid under
          arbitrary heteroskedasticity for errors that are independent across
          cells and symmetric about zero, but NOT under additive cluster
          effects eta_i + xi_j (see ?mwperm_dyadic_het)
      

---

    Code
      print(mwperm_check(index = pn[c("i", "j", "t")]))
    Output
      
      mwperm design diagnosis
      ------------------------------
      Detected design : panel ('t' identified as time by name)
      Roles           : row = i, col = j, time = t
      Dimensions      : i (20) x j (20) x t (3) | 1200 observations
      Balance         : complete
      Resolution      : default K = 19, so p-values are multiples of 1/20 = 0.050
                        -> fine enough for a 95% confidence set at alpha = 0.05
      Would run       : mwperm_panel(y, d, x, row = i, col = j, time = t, time_fe = TRUE)
      

---

    Code
      print(mwperm_check(index = pm[c("i", "j", "t")]))
    Output
      
      mwperm design diagnosis
      ------------------------------
      Detected design : panel_missing ('t' identified as time by name; incomplete array)
      Roles           : row = i, col = j, time = t
      Dimensions      : i (20) x j (20) x t (3) | 1170 observations
      Balance         : incomplete
      Resolution      : set by the biclique blocks, so not known until
                        they are found (see find_bicliques)
      Would run       : mwperm_panel_missing(y, d, x, row = i, col = j, time = t, L0 = NULL, min_block = ..., time_fe = TRUE)
        - The array is incomplete, so the test restricts to (row, col) pairs
          observed in EVERY period and to the fully observed blocks the
          biclique search extracts from them; the group order follows those
          blocks. If some pairs miss periods, `L0 =` keeps instead the L0
          periods jointly observed by the most pairs, the same ones in every
          cell. See ?mwperm_panel_missing and find_bicliques().
      

---

    Code
      print(mwperm_check(index = ly[c("i", "j")], rep = ly$l))
    Output
      
      mwperm design diagnosis
      ------------------------------
      Detected design : layout (`rep =` declares within-cell replication)
      Roles           : row = i, col = j, rep = rep
      Dimensions      : i (4) x j (4) | 113 observations
      Balance         : replicated (16 cells, 6-8 replicates; the smallest cell sets K)
      Resolution      : default K = 5, so p-values are multiples of 1/6 = 0.167
                        -> TOO COARSE for a 95% confidence set at alpha = 0.05 (p
                           cannot reach 0.05). The p-value is still exact; a 95% set
                           needs >= 20 levels in the smallest permuted dimension.
      Would run       : mwperm_layout(y, d, x, row = i, col = j, rep = rep)
      

---

    Code
      print(mwperm_check(index = tw[c("i", "j", "l")], design = "threeway"))
    Output
      
      mwperm design diagnosis
      ------------------------------
      Detected design : threeway (forced via design =)
      Roles           : id1 = i, id2 = j, id3 = l
      Dimensions      : i (6) x j (6) x l (6) | 216 observations
      Balance         : complete
      Resolution      : default K = 5, so p-values are multiples of 1/6 = 0.167
                        -> TOO COARSE for a 95% confidence set at alpha = 0.05 (p
                           cannot reach 0.05). The p-value is still exact; a 95% set
                           needs >= 20 levels in the smallest permuted dimension.
      Would run       : mwperm_threeway(y, d, x, id1 = i, id2 = j, id3 = l)
      

---

    Code
      print(mwperm_check(index = ly[c("i", "j")], design = "irregular"))
    Output
      
      mwperm design diagnosis
      ------------------------------
      Detected design : irregular (forced via design =)
      Roles           : row = i, col = j, rep = (within-cell order)
      Dimensions      : i (4) x j (4) | 113 observations
      Balance         : irregular (16 cells, 6-8 observations each; L0 sets which cells are usable)
      Resolution      : set by the biclique blocks, so not known until
                        they are found (see find_bicliques)
      Would run       : mwperm_irregular(y, d, x, row = i, col = j, L0 = ...)
        - The permutation-group order for the Section 6.4 design is set by
          the biclique blocks found under the mask M_ij = 1{ell_ij >= L0}, so
          it depends on `L0`; see find_bicliques() and ?mwperm_irregular. Its
          random per-cell trim is exact only for exchangeable replicates
          within a cell; if the within-cell index is a period, or anything
          else with an effect shared across cells, it can over-reject even
          with a cell-constant `d` -- use `time =` with `L0 =` instead
          (mwperm_panel_missing(), the same L0 periods in every cell).
      

---

    Code
      print(mwperm_check(index = dy[c("i", "j")], design = "dyadic_het"))
    Output
      
      mwperm design diagnosis
      ------------------------------
      Detected design : dyadic_het (forced via design =)
      Roles           : row = i, col = j
      Dimensions      : i (20) x j (20) | 400 observations
      Balance         : complete
      Resolution      : default n_flip = 6, so p-values are multiples of 1/2^5 = 1/32 = 0.031
                        -> fine enough for a 95% confidence set at alpha = 0.05
      Would run       : mwperm_dyadic_het(y, d, x, row = i, col = j)
      

