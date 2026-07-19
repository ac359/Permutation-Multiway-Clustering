## Phase 2.1 -- build_perm_set() is Algorithm 1 (GTW 2026): a random
## block-cyclic
## permutation group. These are EXACT invariants (machine precision /
## identical),
## the class of bug Monte Carlo will not reliably surface. Base-R stopifnot so
## it
## ships and runs under R CMD check without a testthat dependency.
library(mwperm)

## composition of image vectors: (a o b)(i) = a[b[i]]
compose <- function(a, b) a[b]
## unpack a c(n=, K=) config as integers (a bare c() would coerce to double and
## break identical() against the integer block_size / group indices)
cfg_nK <- function(cfg)
  list(n = as.integer(cfg[["n"]]), K = as.integer(cfg[["K"]]))

## ---- 1. each element is a bijection of 1:n; element 1 is the identity -------
for (cfg in list(list(n = 12, K = 3), list(n = 8, K = 7), list(n = 25,
                                                               K = 4))) {
  p <- cfg_nK(cfg)
  n <- p$n
  K <- p$K
  G <- build_perm_set(n, K, seed = 1)
  stopifnot(length(G) == K + 1L)
  stopifnot(identical(G[[1]], seq_len(n)))                    # identity first
  for (g in G) stopifnot(identical(sort(g), seq_len(n)))      # bijection of 1:n
  stopifnot(identical(attr(G, "block_size"), K + 1L))
}

## ---- 2. closure: the WHOLE point (Proposition 2 / brief H4) -----------------
## Build the full (K+1)x(K+1) composition table; every composite must be the
## group element indexed by (r + s) mod (K+1). Equivalently g_k = g_1^k (cyclic,
## one generator) -- the property Theorem 1's proof relies on.
for (cfg in list(list(n = 12, K = 3), list(n = 20, K = 4), list(n = 9,
                                                                K = 2))) {
  p <- cfg_nK(cfg)
  n <- p$n
  K <- p$K
  B <- K + 1L
  G <- build_perm_set(n, K, seed = 7)
  for (r in 0:K) for (s in 0:K)
    stopifnot(identical(compose(G[[r + 1L]], G[[s + 1L]]),
                        G[[((r + s) %% B) + 1L]]))            # closure
  it <- seq_len(n)                                       # g_1^k accumulator
  for (k in 1:K) {
    it <- compose(G[[2L]], it)
    stopifnot(identical(it, G[[k + 1L]]))
  }
}

## ---- 3. remainder handling (brief H4): leftover tail is fixed, no interior FP
## -
## When n %% (K+1) != 0 the leftover indices cannot be cyclically shifted within
## a
## short block without breaking closure, so they MUST be fixed points; every
## in-block index MUST move under every non-identity element.
for (cfg in list(list(n = 10, K = 3), list(n = 7, K = 3), list(n = 25, K = 4),
                 list(n = 13, K = 12), list(n = 6, K = 1))) {
  p <- cfg_nK(cfg)
  n <- p$n
  K <- p$K
  B <- K + 1L
  tail_size <- n %% B                            # expected # fixed points
  G <- build_perm_set(n, K, seed = 3)
  fp_sets <- lapply(2:(K + 1L), function(k) which(G[[k]] == seq_len(n)))
  ## every non-identity element fixes exactly the tail (same set, size n %% B)
  for (fp in fp_sets) stopifnot(length(fp) == tail_size)
  stopifnot(length(unique(fp_sets)) == 1L)          # identical fixed set
}

## ---- 4. input validation: K+1 <= n, K >= 1, n >= 2 --------------------------
err <- function(e)
  tryCatch({
    e
    NA_character_
  }, error = function(x) conditionMessage(x))
stopifnot(!is.na(err(build_perm_set(5, 5))))     # K + 1 = 6 > n = 5
stopifnot(!is.na(err(build_perm_set(5, 0))))     # K < 1
stopifnot(!is.na(err(build_perm_set(1, 1))))     # n < 2

## ---- 5. RNG hygiene: a seeded call leaves .Random.seed byte-identical -------
## (brief 2.1) test under both RNG kinds, and the case where .Random.seed does
## not yet exist.
for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
  RNGkind(kind)
  set.seed(999)
  before <- .Random.seed
  invisible(build_perm_set(30, 5, seed = 1))
  stopifnot(identical(.Random.seed, before))     # caller stream untouched
}
RNGkind("Mersenne-Twister")
## no pre-existing .Random.seed: a seeded call must not leave one behind
if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed",
                                                   envir = .GlobalEnv)
invisible(build_perm_set(30, 5, seed = 1))
stopifnot(!exists(".Random.seed", envir = .GlobalEnv))

## ---- 6. determinism: same seed -> same group; distinct seeds -> distinct ----
stopifnot(identical(build_perm_set(40, 6, seed = 2), build_perm_set(40, 6,
                                                                    seed = 2)))
stopifnot(!identical(build_perm_set(40, 6, seed = 2), build_perm_set(40, 6,
                                                                     seed = 3)))

cat("test-permset.R: all assertions passed\n")
