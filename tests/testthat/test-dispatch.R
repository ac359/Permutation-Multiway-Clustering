## The unified entry point: every row of the JSS draft's Table 1 (Section
## 3.2) routes to the worker it names, mwperm() returns what that worker
## returns for the same seed, a complete three-index array defaults to the
## panel design (three-way only on request), and mwperm_check() prints the
## diagnosis for each design (snapshots).

same_numbers <- function(a, b) {
  expect_identical(a$pvalues_rep, b$pvalues_rep)
  expect_identical(a$pvalue, b$pvalue)
  expect_identical(unname(a$estimate), unname(b$estimate))
  expect_identical(a$conf_int, b$conf_int)
  expect_identical(a$conf_set, b$conf_set)
  expect_identical(a$K, b$K)
  expect_identical(a$n_obs, b$n_obs)
}
quiet <- function(expr) suppressMessages(suppressWarnings(expr))

dy <- make_dyadic(20, 20, beta = 0.2, seed = 101)
md <- dy[with_seed(102, stats::runif(nrow(dy)) < 0.85), ]
pn <- make_panel(20, 20, T = 3, beta = 0.2, seed = 103)
pm <- pn[!(paste(pn$i, pn$j) %in% with_seed(104, sample(
  unique(paste(pn$i, pn$j)), 15)) & pn$t > 1), ]
ly <- make_layout(4, 4, sizes = 6:8, beta = 0.3, seed = 105)
tw <- make_threeway(6, 6, 6, beta = 0.2, seed = 106)

test_that("Table 1: two indices, complete array -> mwperm_dyadic()", {
  expect_identical(mwperm_check(index = dy[c("i", "j")])$design, "dyadic")
  a <- quiet(mwperm(dy$y, dy$d, x = dy$x, index = dy[c("i", "j")],
                    n_reps = 3, seed = 1))
  expect_identical(a$auto$design, "dyadic")
  same_numbers(a, mwperm_dyadic(dy$y, dy$d, x = dy$x, row = dy$i,
                                col = dy$j, n_reps = 3, seed = 1))
})

test_that("Table 1: two indices, cells missing -> mwperm_missing()", {
  expect_identical(mwperm_check(index = md[c("i", "j")])$design, "missing")
  a <- quiet(mwperm(md$y, md$d, x = md$x, index = md[c("i", "j")],
                    min_block = 3, n_reps = 3, seed = 1, conf_int = FALSE))
  expect_identical(a$auto$design, "missing")
  same_numbers(a, mwperm_missing(md$y, md$d, x = md$x, row = md$i,
                                 col = md$j, min_block = 3, n_reps = 3,
                                 seed = 1, conf_int = FALSE))
})

test_that("Table 1: three indices (or time =), complete -> mwperm_panel()", {
  expect_identical(mwperm_check(index = pn[c("i", "j", "t")])$design,
                   "panel")
  expect_identical(mwperm_check(index = pn[c("i", "j")], time = pn$t)$design,
                   "panel")
  direct <- mwperm_panel(pn$y, pn$d, x = pn$x, row = pn$i, col = pn$j,
                         time = pn$t, n_reps = 3, seed = 1)
  a <- quiet(mwperm(pn$y, pn$d, x = pn$x, index = pn[c("i", "j", "t")],
                    n_reps = 3, seed = 1))
  b <- quiet(mwperm(pn$y, pn$d, x = pn$x, index = pn[c("i", "j")],
                    time = pn$t, n_reps = 3, seed = 1))
  expect_identical(a$auto$design, "panel")
  same_numbers(a, direct)
  same_numbers(b, direct)
})

test_that("Table 1: three indices (or time =), incomplete -> mwperm_panel_missing()", {
  expect_identical(mwperm_check(index = pm[c("i", "j", "t")])$design,
                   "panel_missing")
  a <- quiet(mwperm(pm$y, pm$d, x = pm$x, index = pm[c("i", "j")],
                    time = pm$t, min_block = 3, n_reps = 3, seed = 1,
                    conf_int = FALSE))
  expect_identical(a$auto$design, "panel_missing")
  same_numbers(a, mwperm_panel_missing(pm$y, pm$d, x = pm$x, row = pm$i,
                                       col = pm$j, time = pm$t,
                                       min_block = 3, n_reps = 3, seed = 1,
                                       conf_int = FALSE))
})

test_that("Table 1: two indices plus rep -> mwperm_layout()", {
  expect_identical(mwperm_check(index = ly[c("i", "j")], rep = ly$l)$design,
                   "layout")
  a <- quiet(mwperm(ly$y, ly$d, x = ly$x, index = ly[c("i", "j")],
                    rep = ly$l, n_reps = 3, seed = 1, conf_int = FALSE))
  expect_identical(a$auto$design, "layout")
  same_numbers(a, mwperm_layout(ly$y, ly$d, x = ly$x, row = ly$i,
                                col = ly$j, rep = ly$l, n_reps = 3, seed = 1,
                                conf_int = FALSE))
})

test_that("Table 1: design = 'threeway' -> mwperm_threeway(), only on request", {
  ## a complete three-index array whose indices all look exchangeable still
  ## defaults to the panel design, which stays valid if one of them is time
  expect_identical(quiet(mwperm_check(index = tw[c("i", "j", "l")]))$design,
                   "panel")
  expect_identical(mwperm_check(index = tw[c("i", "j", "l")],
                                design = "threeway")$design, "threeway")
  a <- quiet(mwperm(tw$y, tw$d, x = tw$x, index = tw[c("i", "j", "l")],
                    design = "threeway", n_reps = 3, seed = 1,
                    conf_int = FALSE))
  expect_identical(a$auto$design, "threeway")
  same_numbers(a, mwperm_threeway(tw$y, tw$d, x = tw$x, id1 = tw$i,
                                  id2 = tw$j, id3 = tw$l, n_reps = 3,
                                  seed = 1, conf_int = FALSE))
})

test_that("Table 1: design = 'irregular' with L0 -> mwperm_irregular()", {
  li <- make_layout(8, 8, sizes = 3:6, seed = 107, d_cell_constant = TRUE)
  a <- quiet(mwperm(li$y, li$d, x = li$x, index = li[c("i", "j")],
                    rep = li$l, design = "irregular", L0 = 4, min_block = 3,
                    n_reps = 3, seed = 1, conf_int = FALSE))
  expect_identical(a$auto$design, "irregular")
  same_numbers(a, mwperm_irregular(li$y, li$d, x = li$x, row = li$i,
                                   col = li$j, rep = li$l, L0 = 4,
                                   min_block = 3, n_reps = 3, seed = 1,
                                   conf_int = FALSE))
  ## never inferred, and L0 is required
  expect_identical(quiet(mwperm_check(index = li[c("i", "j")]))$design,
                   "layout")
  expect_error(quiet(mwperm(li$y, li$d, index = li[c("i", "j")],
                            design = "irregular", seed = 1)), "requires `L0")
})

test_that("Table 1: design = 'dyadic_het' -> mwperm_dyadic_het(), only on request", {
  expect_identical(mwperm_check(index = dy[c("i", "j")],
                                design = "dyadic_het")$design, "dyadic_het")
  a <- quiet(mwperm(dy$y, dy$d, x = dy$x, index = dy[c("i", "j")],
                    design = "dyadic_het", n_flip = 5, n_reps = 3, seed = 1,
                    conf_int = FALSE))
  expect_identical(a$auto$design, "dyadic_het")
  same_numbers(a, mwperm_dyadic_het(dy$y, dy$d, x = dy$x, row = dy$i,
                                    col = dy$j, n_flip = 5, n_reps = 3,
                                    seed = 1, conf_int = FALSE))
  ## K means nothing to a sign-flip group: warned and ignored
  expect_warning(suppressMessages(mwperm(
    dy$y, dy$d, index = dy[c("i", "j")], design = "dyadic_het", K = 5,
    n_reps = 1, seed = 1, conf_int = FALSE)), "does not apply")
})

test_that("mwperm_check() prints a diagnosis for every design", {
  local_reproducible_output(width = 80)
  expect_snapshot(print(mwperm_check(index = dy[c("i", "j")])))
  expect_snapshot(print(mwperm_check(index = md[c("i", "j")])))
  expect_snapshot(print(mwperm_check(index = pn[c("i", "j", "t")])))
  expect_snapshot(print(mwperm_check(index = pm[c("i", "j", "t")])))
  expect_snapshot(print(mwperm_check(index = ly[c("i", "j")], rep = ly$l)))
  expect_snapshot(print(mwperm_check(index = tw[c("i", "j", "l")],
                                     design = "threeway")))
  expect_snapshot(print(mwperm_check(index = ly[c("i", "j")],
                                     design = "irregular")))
  expect_snapshot(print(mwperm_check(index = dy[c("i", "j")],
                                     design = "dyadic_het")))
})
