## Tests for build_perm_set(): Algorithm 1 produces a valid cyclic group.
## Plain base-R tests (testthat is not assumed available).

library(mwperm)

## helper: compose two image-vector permutations (apply p1, then p2)
compose <- function(p2, p1) p2[p1]

set.seed(42)

for (trial in seq_len(5)) {
  n <- sample(6:20, 1)
  K <- sample(2:(n - 1), 1)
  G <- build_perm_set(n, K, seed = trial)

  ## (a) correct count and identity first element
  stopifnot(length(G) == K + 1L)
  stopifnot(identical(G[[1]], seq_len(n)))
  stopifnot(attr(G, "block_size") == K + 1L)

  ## (b) every element is a genuine permutation of seq_len(n)
  for (g in G) stopifnot(sort(g) == seq_len(n))

  ## (c) closure: composing any two elements lands back in the group
  keys <- vapply(G, paste, character(1), collapse = ",")
  for (a in seq_along(G)) for (b in seq_along(G)) {
    comp <- compose(G[[b]], G[[a]])
    stopifnot(paste(comp, collapse = ",") %in% keys)
  }

  ## (d) the non-identity generator has order K+1 (cyclic group)
  gen <- G[[2]]
  pw <- seq_len(n)
  for (e in seq_len(K)) {
    pw <- compose(gen, pw)
    stopifnot(!identical(pw, seq_len(n)))   # not identity before order K+1
  }
  pw <- compose(gen, pw)
  stopifnot(identical(pw, seq_len(n)))        # identity at power K+1
}

## error handling: K+1 must not exceed n
ok <- tryCatch({ build_perm_set(4, 4); FALSE }, error = function(e) TRUE)
stopifnot(ok)

cat("test-perm-set.R: all checks passed\n")
