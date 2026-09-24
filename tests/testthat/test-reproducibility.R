## Reproducibility, as the JSS draft states it (Section 3.8): a seed fixes the
## result completely; repetition r draws from seed + r - 1; parallel runs are
## bit-for-bit identical to serial ones; and a seeded fit leaves the caller's
## random stream where it found it. Whether an UNSEEDED fit moves the global
## stream is recorded as well (it does, by design: it draws from it).

numeric_fields <- c("pvalue", "pvalues_rep", "estimate", "se_naive",
                    "conf_int", "conf_set", "conf_region", "K", "n_obs")
same_fit <- function(a, b)
  for (nm in numeric_fields) expect_identical(a[[nm]], b[[nm]], info = nm)

dy <- make_dyadic(20, 20, beta = 0.2, seed = 131)
pn <- make_panel(20, 20, T = 3, beta = 0.2, seed = 132)
ly <- make_layout(4, 4, sizes = 6:9, beta = 0.3, seed = 133)
li <- make_layout(8, 8, sizes = 3:6, seed = 134, d_cell_constant = TRUE)
md <- subset(dy, i != j)
tw <- make_threeway(6, 6, 7, seed = 135)

fitters <- list(
  dyadic = function(...) mwperm_dyadic(dy$y, dy$d, x = dy$x, row = dy$i,
                                       col = dy$j, ...),
  panel = function(...) mwperm_panel(pn$y, pn$d, x = pn$x, row = pn$i,
                                     col = pn$j, time = pn$t, ...),
  threeway = function(...) mwperm_threeway(tw$y, tw$d, id1 = tw$i,
                                           id2 = tw$j, id3 = tw$l, ...),
  layout_L0 = function(...) mwperm_layout(ly$y, ly$d, row = ly$i,
                                          col = ly$j, rep = ly$l, L0 = 7,
                                          ...),
  missing = function(...) mwperm_missing(md$y, md$d, row = md$i,
                                         col = md$j, ...),
  irregular = function(...) mwperm_irregular(li$y, li$d, row = li$i,
                                             col = li$j, rep = li$l, L0 = 4,
                                             min_block = 3, ...),
  sign_flip = function(...) mwperm_dyadic_het(dy$y, dy$d, row = dy$i,
                                              col = dy$j, n_flip = 5, ...))

test_that("the same seed gives an identical fit, in every design", {
  for (nm in names(fitters)) {
    f <- fitters[[nm]]
    same_fit(f(n_reps = 3, seed = 17), f(n_reps = 3, seed = 17))
  }
})

test_that("repetition r draws from seed + r - 1", {
  for (nm in c("dyadic", "panel", "missing", "sign_flip")) {
    f <- fitters[[nm]]
    all3 <- f(n_reps = 3, seed = 40, conf_int = FALSE)$pvalues_rep
    one_at_a_time <- vapply(40:42, function(s)
      f(n_reps = 1, seed = s, conf_int = FALSE)$pvalues_rep, 0)
    expect_identical(all3, one_at_a_time, info = nm)
  }
})

test_that("n_cores = 2 is identical() to n_cores = 1", {
  skip_on_cran()
  skip_on_os("windows")   # PSOCK workers need an installed package
  skip_if(parallel::detectCores() < 2L, "fewer than two cores")
  for (nm in names(fitters)) {
    f <- fitters[[nm]]
    ## repetition axis (seeded, n_reps > 1) and the per-element axis
    same_fit(f(n_reps = 4, seed = 9, n_cores = 2),
             f(n_reps = 4, seed = 9, n_cores = 1))
    same_fit(f(n_reps = 1, seed = 9, n_cores = 2),
             f(n_reps = 1, seed = 9, n_cores = 1))
  }
})

test_that("a seeded fit leaves the global RNG state untouched", {
  set.seed(2026)
  for (nm in names(fitters)) {
    before <- .Random.seed
    fitters[[nm]](n_reps = 2, seed = 3)
    expect_identical(.Random.seed, before, info = nm)
  }
  ## also when no RNG state exists yet
  if (exists(".Random.seed", envir = globalenv()))
    rm(".Random.seed", envir = globalenv())
  fitters$dyadic(n_reps = 1, seed = 3, conf_int = FALSE)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  set.seed(1)
})

test_that("an unseeded fit draws from, and so advances, the global stream", {
  ## Recorded behaviour: seed = NULL means the caller's stream IS the source
  ## of the relabellings, so it moves; re-seeding reproduces the fit.
  set.seed(77)
  before <- .Random.seed
  a <- fitters$dyadic(n_reps = 2, conf_int = FALSE)
  expect_false(identical(.Random.seed, before))
  set.seed(77)
  b <- fitters$dyadic(n_reps = 2, conf_int = FALSE)
  expect_identical(a$pvalues_rep, b$pvalues_rep)
})
