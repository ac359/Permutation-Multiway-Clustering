## mwperm() and mwperm_check() -- the unified entry point and the detector.
##
## mwperm() is a dispatcher: it classifies the clustering structure, announces
## what it found, and calls the matching design-specific front end. Two things
## must therefore hold. Dispatch has to be an IDENTITY -- the same seed through
## mwperm() and through the front end directly must give the same object, field
## for field -- and every fork the detector takes on incomplete evidence has to
## be announced, because the wrong design is the largest practical risk with a
## finite-sample-valid test: treating autocorrelated time as an exchangeable
## third dimension is a silent size distortion, not an error message.
##
## The per-design behaviour is in test-dyadic.R, test-panel.R and friends; see
## tests/README.md.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

## ---- fixtures -------------------------------------------------------------
set.seed(202)
g2 <- expand.grid(i = 1:6, j = 1:6)                 # complete two-way
y2 <- rnorm(36)
d2 <- rnorm(36)
df2 <- data.frame(i = g2$i, j = g2$j, yy = y2, dd = d2, xx = rnorm(36))
g3 <- expand.grid(i = 1:5, j = 1:5, year = 1:3)     # three indices
gm <- g2[g2$i != g2$j, ]                            # incomplete two-way

## ---- 1. what the detector reports, and how it prints ----------------------
chk <- mwperm_check(index = list(i = g2$i, j = g2$j))
stopifnot(inherits(chk, "mwperm_design"), identical(chk$design, "dyadic"))
out_d <- paste(capture.output(print(chk)), collapse = "\n")
stopifnot(grepl("dyadic", out_d), grepl("Resolution", out_d),
          grepl("TOO COARSE", out_d))               # 6 clusters, needs 20
out_m <- paste(capture.output(print(
  mwperm_check(index = list(i = gm$i, j = gm$j)))), collapse = "\n")
stopifnot(grepl("missing", out_m), grepl("biclique", out_m))
out_p <- paste(capture.output(print(
  mwperm_check(index = list(i = g3$i, j = g3$j, year = g3$year)))),
  collapse = "\n")
stopifnot(grepl("panel", out_p))
## a diagnosis carrying warnings renders them with the "!" prefix
out_w <- paste(capture.output(print(
  mwperm_check(index = list(period = g3$i, a = g3$j, b = g3$year)))),
  collapse = "\n")
stopifnot(grepl("!", out_w, fixed = TRUE))
## forced missing on an incomplete array: K depends on blocks not yet searched
chk_fm <- mwperm_check(index = list(i = gm$i, j = gm$j), design = "missing")
stopifnot(identical(chk_fm$design, "missing"), is.na(chk_fm$K_default))

## The resolution verdict is about the alpha and aggregation the fit will use,
## not a hard-coded 0.05 under the default rule. K = 5 here: the floor 1/6 is
## too coarse at 0.05 but fine at 0.2 -- and under "median2", whose floor is
## 2/6, too coarse at 0.2 as well. The printed line names the alpha.
chk_a <- mwperm_check(index = list(i = g2$i, j = g2$j), alpha = 0.2)
chk_m <- mwperm_check(index = list(i = g2$i, j = g2$j), alpha = 0.2,
                      aggregate = "median2")
stopifnot(identical(chk$resolution_ok, FALSE), identical(chk$alpha, 0.05),
          identical(chk_a$resolution_ok, TRUE), identical(chk_a$alpha, 0.2),
          identical(chk_m$resolution_ok, FALSE),
          identical(chk$p_floor, 1 / 6), identical(chk_m$p_floor, 2 / 6),
          chk$levels_needed == 20L, chk_a$levels_needed == 5L,
          chk_m$levels_needed == 10L)
out_a <- paste(capture.output(print(chk_a)), collapse = "\n")
out_m2 <- paste(capture.output(print(chk_m)), collapse = "\n")
stopifnot(grepl("alpha = 0.2", out_a, fixed = TRUE),
          grepl("80% confidence set", out_a, fixed = TRUE),
          !grepl("TOO COARSE", out_a),
          grepl("TOO COARSE", out_m2), grepl("median2", out_m2, fixed = TRUE),
          grepl("alpha = 0.05", out_d, fixed = TRUE))
expect_err(mwperm_check(index = list(i = g2$i, j = g2$j), alpha = 1),
           "`alpha`")

## ---- 2. index resolution and forced-design validation ---------------------
## Forcing a design that the data cannot support must fail with a message about
## the DATA's shape, not a downstream indexing error.
expect_err(mwperm_check(index = c("i", "j")), "`data`")
expect_err(mwperm_check(index = c("i", "nope"), data = df2), "not found")
expect_err(mwperm_check(index = 1:5), "`index`")
expect_err(mwperm_check(index = list(a = g3$i, b = g3$j, cc = g3$year,
                                     dd = g3$year)),
           "Too many clustering dimensions")
expect_err(mwperm_check(index = list(i = g3$i, j = g3$j, t = g3$year),
                        design = "dyadic"), "exactly 2")
dup2 <- rbind(g2, g2)
expect_err(mwperm_check(index = list(i = dup2$i, j = dup2$j),
                        design = "dyadic"), "one observation per cell")
expect_err(mwperm_check(index = list(i = gm$i, j = gm$j), design = "dyadic"),
           "missing")
expect_err(mwperm_check(index = list(i = g3$i, j = g3$j, t = g3$year),
                        design = "missing"), "exactly 2")
expect_err(mwperm_check(index = list(i = g2$i, j = g2$j), design = "panel"),
           "time dimension")
expect_err(mwperm_check(index = list(i = g2$i, j = g2$j), design = "threeway"),
           "exactly 3")
gm3 <- g3[-1, ]                                     # incomplete 3-index array
expect_err(mwperm_check(index = list(i = gm3$i, j = gm3$j, t = gm3$year),
                        design = "threeway"), "complete balanced")
expect_err(mwperm_check(index = list(i = g3$i, j = g3$j, t = g3$year),
                        design = "layout"), "exactly 2")
expect_err(mwperm_check(index = list(i = g3$i, j = g3$j, t = g3$year),
                        time = g3$year), "Too many clustering dimensions")
expect_err(mwperm_check(index = list(i = g3$i, j = g3$j, t = g3$year),
                        rep = g3$year), "Too many clustering dimensions")
## An incomplete array with a tagged `time =` is an incomplete PANEL: it routes
## to mwperm_panel_missing(), exactly as three untagged indices already do,
## and forcing design = "panel_missing" -- with the time role tagged or taken
## from the third index -- lands in the same place with the roles filled in.
## Before 0.4.0's dispatcher caught up, the tagged path errored and the forced
## path returned a diagnosis with no roles and an empty "Would run" line.
chk_t <- mwperm_check(index = list(i = gm3$i, j = gm3$j), time = gm3$year)
chk_f <- mwperm_check(index = list(i = gm3$i, j = gm3$j, year = gm3$year),
                      design = "panel_missing")
chk_ft <- mwperm_check(index = list(i = gm3$i, j = gm3$j), time = gm3$year,
                       design = "panel_missing")
for (ck in list(chk_t, chk_f, chk_ft))
  stopifnot(identical(ck$design, "panel_missing"),
            identical(names(ck$roles), c("row", "col", "time")),
            identical(ck$roles$row, "i"), identical(ck$roles$col, "j"),
            is.na(ck$K_default), identical(ck$balance, "incomplete"),
            grepl("^mwperm_panel_missing\\(y, d, x, row = i, col = j, time = ",
                  ck$call_str))
stopifnot(identical(chk_t$roles$time, "time"),
          identical(chk_f$roles$time, "year"),
          grepl("incomplete array", chk_t$reason, fixed = TRUE),
          any(grepl("observed in EVERY period", chk_t$notes, fixed = TRUE)))
out_ft <- paste(capture.output(print(chk_ft)), collapse = "\n")
stopifnot(grepl("Roles           : row = i, col = j, time = time", out_ft,
                fixed = TRUE),
          grepl("Would run       : mwperm_panel_missing(y, d, x, row = i, ",
                out_ft, fixed = TRUE))
## forcing the complete-array panel on an incomplete array names the way out
expect_err(mwperm_check(index = list(i = gm3$i, j = gm3$j), time = gm3$year,
                        design = "panel"), "panel_missing")
## the incomplete panel still needs a time role to hold fixed
expect_err(mwperm_check(index = list(i = g2$i, j = g2$j),
                        design = "panel_missing"), "time dimension")
## a repeated (row, col, time) cell is fatal under a tagged time, as before
dup3 <- rbind(g3, g3[1L, ])
expect_err(mwperm_check(index = list(i = dup3$i, j = dup3$j),
                        time = dup3$year), "repeat")

## ---- 3. detection safety: every uncertain fork is announced ---------------
## A complete balanced 6 x 6 x 4 array can be read as three-way (permute all
## three) or as a panel (hold one fixed). The two differ in what they assume,
## so whenever the evidence for the choice is thin the user has to hear about
## it -- and, just as importantly, must NOT be warned when the evidence is
## clear, or the warnings stop being read.
set.seed(42)
gx <- expand.grid(a = 1:6, b = 1:6, tt = 1:4)
Nx <- nrow(gx)
yx <- rnorm(Nx)
dx <- rnorm(Nx)

## 3a. name-alone evidence: 'period' is a genuine 6-level cluster while 'g' has
## strictly the fewest levels. The name still wins the role -- dispatch is
## frozen -- but it must carry a warning with override instructions.
idx_name <- list(period = gx$a, i = gx$b, g = gx$tt)
c1 <- mwperm_check(index = idx_name)
stopifnot(identical(c1$design, "panel"), identical(c1$roles$time, "period"),
          any(grepl("name", c1$warnings, ignore.case = TRUE) &
              grepl("time =", c1$warnings, fixed = TRUE)))
w1 <- warns_of(f1 <- mwperm(y = yx, d = dx, index = idx_name,
                            conf_int = FALSE, seed = 1, verbose = FALSE))
stopifnot(inherits(f1, "mwperm"),
          any(grepl("name", w1, ignore.case = TRUE)),
          any(grepl("name", f1$note, ignore.case = TRUE)))   # kept on the fit

## 3b. corroborated evidence stays silent: 'year' is regularly spaced with
## strictly fewer levels than both clusters, so the values single it out
idx_ok <- list(i = gx$a, j = gx$b, year = 2000L + gx$tt)
c2 <- mwperm_check(index = idx_ok)
stopifnot(identical(c2$design, "panel"), identical(c2$roles$time, "year"),
          length(c2$warnings) == 0L,
          length(warns_of(mwperm(y = yx, d = dx, index = idx_ok,
                                 conf_int = FALSE, seed = 1,
                                 verbose = FALSE))) == 0L)

## 3c. forcing threeway over a time-like index warns; forcing it where nothing
## looks time-like does not
c3 <- mwperm_check(index = idx_ok, design = "threeway")
stopifnot(identical(c3$design, "threeway"),
          any(grepl("time-like", c3$warnings) & grepl("panel", c3$warnings)))
w3 <- warns_of(f3 <- mwperm(y = yx, d = dx, index = idx_ok,
                            design = "threeway", conf_int = FALSE, seed = 1,
                            verbose = FALSE))
stopifnot(inherits(f3, "mwperm"), any(grepl("time-like", w3)))
g5 <- expand.grid(a = 1:5, b = 1:5, cc = 1:5)       # equal levels: no evidence
set.seed(7)
y5 <- rnorm(nrow(g5))
d5 <- rnorm(nrow(g5))
idx_ex <- list(a = g5$a, b = g5$b, cc = g5$cc)
stopifnot(length(mwperm_check(index = idx_ex, design = "threeway")$warnings)
            == 0L,
          length(warns_of(mwperm(y = y5, d = d5, index = idx_ex,
                                 design = "threeway", conf_int = FALSE,
                                 seed = 1, verbose = FALSE))) == 0L)

## 3d. no evidence either way: the panel default is taken and the warning must
## say the default protects only if the PERMUTED pair is exchangeable -- not
## that it is valid whatever the three dimensions are
c5 <- mwperm_check(index = list(a = gx$a, b = gx$b, cc = gx$tt))
stopifnot(identical(c5$design, "panel"), identical(c5$roles$time, "cc"),
          any(grepl("permuted", c5$warnings) &
              grepl("time =", c5$warnings, fixed = TRUE)))

## 3e. replicated cells: the layout fork is announced at fit time
gr <- expand.grid(i = 1:5, j = 1:5)
gr <- rbind(gr, gr, gr)                             # three replicates per cell
set.seed(11)
expect_warn(mwperm(y = rnorm(nrow(gr)), d = rnorm(nrow(gr)),
                   index = list(i = gr$i, j = gr$j), conf_int = FALSE,
                   seed = 1, verbose = FALSE),
            "replication")

## ---- 4. dispatch is an identity, for every design -------------------------
## Promised by the vignette, the man pages and NEWS: mwperm() must return the
## same object as the direct call, modulo `call` and the `auto` metadata it
## adds. same_fit() reports the first field that differs, so a failure names it.
set.seed(1)
n <- 21L
ge <- expand.grid(i = seq_len(n), j = seq_len(n))
Ne <- nrow(ge)
xe <- cbind(z = rnorm(Ne))
de <- rnorm(n)[ge$i] + rnorm(Ne)
ye <- rnorm(n)[ge$i] + rnorm(n)[ge$j] + 0.4 * de + rnorm(Ne)
dy_dir <- mwperm_dyadic(ye, de, x = xe, row = ge$i, col = ge$j, seed = 5)
dy_dis <- mwperm(y = ye, d = de, x = xe, index = list(row = ge$i, col = ge$j),
                 design = "dyadic", seed = 5, verbose = FALSE)
stopifnot(isTRUE(same_fit(dy_dir, dy_dis, skip = c("call", "auto"))),
          identical(dy_dis$auto$design, "dyadic"))

gp <- expand.grid(i = 1:6, j = 1:6, t = 1:4)
Np <- nrow(gp)
set.seed(2)
dp <- rnorm(Np)
yp <- rnorm(6)[gp$i] + rnorm(6)[gp$j] + cumsum(rnorm(4))[gp$t] + 0.3 * dp +
  rnorm(Np)
pa_dir <- mwperm_panel(yp, dp, row = gp$i, col = gp$j, time = gp$t, seed = 3)
pa_dis <- mwperm(y = yp, d = dp, index = list(row = gp$i, col = gp$j),
                 time = gp$t, design = "panel", seed = 3, verbose = FALSE)
stopifnot(isTRUE(same_fit(pa_dir, pa_dis, skip = c("call", "auto"))))

gt <- expand.grid(a = 1:5, b = 1:5, c = 1:5)
Nt <- nrow(gt)
set.seed(3)
dt3 <- rnorm(Nt)
yt3 <- rnorm(5)[gt$a] + rnorm(5)[gt$b] + rnorm(5)[gt$c] + rnorm(Nt)
tw_dir <- mwperm_threeway(yt3, dt3, id1 = gt$a, id2 = gt$b, id3 = gt$c,
                          seed = 4)
tw_dis <- mwperm(y = yt3, d = dt3,
                 index = list(id1 = gt$a, id2 = gt$b, id3 = gt$c),
                 design = "threeway", seed = 4, verbose = FALSE)
stopifnot(isTRUE(same_fit(tw_dir, tw_dis, skip = c("call", "auto"))))

gl <- expand.grid(l = 1:5, i = 1:4, j = 1:4)
Nl <- nrow(gl)
set.seed(4)
dl <- rnorm(Nl)
yl <- rnorm(16)[as.integer(interaction(gl$i, gl$j))] + 0.5 * dl + rnorm(Nl)
la_dir <- mwperm_layout(yl, dl, row = gl$i, col = gl$j, seed = 6)
la_dis <- suppressWarnings(
  mwperm(y = yl, d = dl, index = list(row = gl$i, col = gl$j),
         design = "layout", seed = 6, verbose = FALSE))
stopifnot(isTRUE(same_fit(la_dir, la_dis, skip = c("call", "auto"))))

gmm <- expand.grid(i = 1:12, j = 1:12)
gmm <- gmm[gmm$i != gmm$j, ]
set.seed(5)
gmm <- gmm[sample(nrow(gmm), round(0.9 * nrow(gmm))), ]
Nm <- nrow(gmm)
dmm <- rnorm(Nm)
ymm <- rnorm(12)[gmm$i] + rnorm(12)[gmm$j] + rnorm(Nm)
mi_dir <- mwperm_missing(ymm, dmm, row = gmm$i, col = gmm$j, min_block = 3,
                         seed = 7)
mi_dis <- mwperm(y = ymm, d = dmm, index = list(row = gmm$i, col = gmm$j),
                 design = "missing", min_block = 3, seed = 7, verbose = FALSE)
stopifnot(isTRUE(same_fit(mi_dir, mi_dis, skip = c("call", "auto"))))

## Incomplete panel, reached three ways: auto-detected from a tagged `time =`,
## auto-detected from three untagged indices, and forced. All must be the
## direct mwperm_panel_missing() call. The auto routes add the detector's note
## to `note`; the forced route adds nothing, so it is compared in full.
gpm <- gp[-c(3L, 40L), ]                            # 6 x 6 x 4 minus 2 cells
set.seed(6)
dpm <- rnorm(nrow(gpm))
ypm <- rnorm(6)[gpm$i] + rnorm(6)[gpm$j] + cumsum(rnorm(4))[gpm$t] +
  0.3 * dpm + rnorm(nrow(gpm))
pm_dir <- mwperm_panel_missing(ypm, dpm, row = gpm$i, col = gpm$j,
                               time = gpm$t, min_block = 3, seed = 8)
pm_tag <- mwperm(y = ypm, d = dpm, index = list(row = gpm$i, col = gpm$j),
                 time = gpm$t, min_block = 3, seed = 8, verbose = FALSE)
pm_3 <- suppressWarnings(
  mwperm(y = ypm, d = dpm, index = list(row = gpm$i, col = gpm$j, t = gpm$t),
         min_block = 3, seed = 8, verbose = FALSE))
pm_frc <- mwperm(y = ypm, d = dpm, index = list(row = gpm$i, col = gpm$j),
                 time = gpm$t, design = "panel_missing", min_block = 3,
                 seed = 8, verbose = FALSE)
stopifnot(isTRUE(same_fit(pm_dir, pm_tag, skip = c("call", "auto", "note"))),
          isTRUE(same_fit(pm_dir, pm_3, skip = c("call", "auto", "note"))),
          isTRUE(same_fit(pm_dir, pm_frc, skip = c("call", "auto"))),
          identical(pm_tag$auto$design, "panel_missing"),
          identical(pm_3$auto$design, "panel_missing"),
          identical(pm_frc$auto$design, "panel_missing"),
          identical(pm_frc$auto$roles,
                    list(row = "row", col = "col", time = "time")),
          any(grepl("observed in EVERY period", pm_tag$note, fixed = TRUE)),
          all(pm_dir$note %in% pm_tag$note))

## `L0` travels with `time =` to the incomplete-panel design (0.4.2): the
## dispatched fit is the direct mwperm_panel_missing(L0 = ) call, it is NOT
## the L0 = NULL fit (a different period set is retained), and no
## "`L0` applies to ... only" warning fires for this design.
gpl <- gp[!(gp$i <= 3L & gp$t == 4L), ]            # rows 1-3 miss period 4
set.seed(7)
dpl <- rnorm(nrow(gpl))
ypl <- rnorm(6)[gpl$i] + rnorm(6)[gpl$j] + c(0, 2, 4, 6)[gpl$t] + 0.3 * dpl +
  rnorm(nrow(gpl))
pl_dir <- mwperm_panel_missing(ypl, dpl, row = gpl$i, col = gpl$j,
                               time = gpl$t, L0 = 3L, min_block = 3,
                               conf_int = FALSE, n_reps = 2L, seed = 8)
wl <- warns_of(pl_dis <- mwperm(y = ypl, d = dpl,
                                index = list(row = gpl$i, col = gpl$j),
                                time = gpl$t, L0 = 3L, min_block = 3,
                                conf_int = FALSE, n_reps = 2L, seed = 8,
                                verbose = FALSE))
pl_null <- mwperm_panel_missing(ypl, dpl, row = gpl$i, col = gpl$j,
                                time = gpl$t, min_block = 3,
                                conf_int = FALSE, n_reps = 2L, seed = 8)
stopifnot(isTRUE(same_fit(pl_dir, pl_dis, skip = c("call", "auto", "note"))),
          identical(pl_dis$auto$design, "panel_missing"),
          identical(pl_dis$periods_used, c("1", "2", "3")),
          pl_dis$cells_used == 36L,                  # every pair has 1-3
          pl_null$cells_used == 18L,                 # only rows 4-6 have 1-4
          !isTRUE(same_fit(pl_dir, pl_null, skip = c("call", "note"))),
          !any(grepl("`L0`", wl, fixed = TRUE)),
          any(grepl("`L0 =`", pl_dis$note, fixed = TRUE)))   # the routing note
chk_l <- mwperm_check(index = list(row = gpl$i, col = gpl$j), time = gpl$t)
stopifnot(grepl("L0 = NULL", chk_l$call_str, fixed = TRUE),
          any(grepl("`L0 =`", chk_l$notes, fixed = TRUE)))

## ---- 5. data = : columns are resolved by name -----------------------------
f_nm <- mwperm(y = "yy", d = c("dd", "xx"), index = c("i", "j"), data = df2,
               seed = 2, n_reps = 1, conf_int = FALSE, verbose = FALSE)
stopifnot(length(f_nm$estimate) == 2L,
          identical(f_nm$d_names, c("dd", "xx")))
expect_err(mwperm(y = "zzz", d = "dd", index = c("i", "j"), data = df2,
                  verbose = FALSE), "`y`")
## verbose = TRUE announces the design it chose and the call it is making
ann <- msgs_of(mwperm(y = "yy", d = "dd", index = c("i", "j"), data = df2,
                      seed = 2, n_reps = 1, conf_int = FALSE, verbose = TRUE))
stopifnot(any(grepl("dyadic", ann)))

## ---- 6. arguments meant for another design warn and are ignored -----------
## Silently dropping them would let a user believe a knob was applied.
w_arg <- warns_of(mwperm(y = "yy", d = "dd", index = c("i", "j"), data = df2,
                         time_fe = FALSE, L0 = 3, seed = 2, n_reps = 1,
                         conf_int = FALSE, verbose = FALSE))
stopifnot(any(grepl("`time_fe`", w_arg, fixed = TRUE)),
          any(grepl("`L0`", w_arg, fixed = TRUE)))
expect_warn(mwperm(y = y2, d = d2, index = list(i = g2$i, j = g2$j),
                   permute = "rows", seed = 1, n_reps = 1, conf_int = FALSE,
                   verbose = FALSE),
            "`permute` applies to the missing design only")
## a list of several designs reads as English, not "A and B and C"
expect_warn(mwperm(y = y2, d = d2, index = list(i = g2$i, j = g2$j),
                   min_block = 3, seed = 1, conf_int = FALSE, verbose = FALSE),
            "applies to the missing, irregular and panel_missing designs only")

## ---- 7. n_reps: the dispatched front end's own default ---------------------
## Since 0.4.2 mwperm_irregular() defaults to 500 repetitions (the paper's
## random per-cell trim is redrawn in every repetition, so the median needs
## many of them) while every other front end keeps 10. mwperm()'s own default
## is therefore NULL, "whatever the front end says", and it forwards n_reps
## only when given: the dispatch identity must hold WITHOUT passing it.
stopifnot(is.null(formals(mwperm)$n_reps),
          identical(formals(mwperm_irregular)$n_reps, 500L),
          identical(formals(mwperm_dyadic)$n_reps, 10L))
set.seed(11)
dat_irr <- do.call(rbind, lapply(1:8, function(i)
  do.call(rbind, lapply(1:8, function(j) {
    L <- sample(c(0L, 2L, 4L, 6L, 9L), 1L)
    if (L == 0L) return(NULL)
    data.frame(i = i, j = j, l = seq_len(L), d = rnorm(1), eta = rnorm(1))
  }))))
dat_irr$y <- 0.5 * dat_irr$d + dat_irr$eta + rnorm(nrow(dat_irr))
ir_dir <- with(dat_irr, mwperm_irregular(y = y, d = d, row = i, col = j,
                                         rep = l, L0 = 4L, min_block = 2L,
                                         conf_int = FALSE, seed = 5))
ir_dis <- mwperm(y = "y", d = "d", index = c("i", "j"), rep = "l",
                 data = dat_irr, design = "irregular", L0 = 4L,
                 min_block = 2L, conf_int = FALSE, seed = 5, verbose = FALSE)
## (the detector prepends its L0 routing note, hence `note` is skipped and
## checked by inclusion, as for the incomplete-panel routes above)
stopifnot(identical(ir_dir$n_reps, 500L), length(ir_dir$pvalues_rep) == 500L,
          identical(ir_dis$n_reps, 500L),
          isTRUE(same_fit(ir_dir, ir_dis, skip = c("call", "auto", "note"))),
          all(ir_dir$note %in% ir_dis$note))
## every other design still runs 10 repetitions through mwperm()
dy10 <- mwperm(y = "yy", d = "dd", index = c("i", "j"), data = df2, seed = 2,
               conf_int = FALSE, verbose = FALSE)
stopifnot(identical(dy10$n_reps, 10L), length(dy10$pvalues_rep) == 10L)
## and an explicit n_reps is forwarded as given, to either
ir3 <- mwperm(y = "y", d = "d", index = c("i", "j"), rep = "l", data = dat_irr,
              design = "irregular", L0 = 4L, min_block = 2L, conf_int = FALSE,
              n_reps = 3L, seed = 5, verbose = FALSE)
stopifnot(identical(ir3$n_reps, 3L))

## ---- 8. design = "irregular" with a `time =` index -> the incomplete panel --
## Section 6.4 with PERIODS as the within-cell index is the incomplete panel
## with one observation per cell and period, and the random trim is exact
## only for exchangeable replicates. A tagged `time =` therefore runs
## mwperm_panel_missing() -- identical to the direct call -- instead of being
## dropped while the trim ran on the order of appearance (the pre-fix
## behaviour). Rows 1-5 observe periods 1-3, rows 6-10 periods 2-4, so no pair
## clears the every-period mask and L0 = 2 keeps periods 2-3 in every cell.
set.seed(12)
st <- do.call(rbind, lapply(1:10, function(i)
  data.frame(i = i, j = rep(1:10, each = 3L),
             t = rep(if (i <= 5L) 1:3 else 2:4, times = 10L))))
st$d <- rnorm(100)[(st$i - 1L) * 10L + st$j]
st$y <- rnorm(10)[st$i] + rnorm(10)[st$j] + 0.5 * st$t + rnorm(nrow(st))
pm_dir <- with(st, mwperm_panel_missing(y = y, d = d, row = i, col = j,
                                        time = t, L0 = 2L, n_reps = 3L,
                                        conf_int = FALSE, seed = 5))
pm_dis <- mwperm(y = "y", d = "d", index = c("i", "j"), time = "t", data = st,
                 design = "irregular", L0 = 2L, n_reps = 3L, conf_int = FALSE,
                 seed = 5, verbose = FALSE)
stopifnot(identical(pm_dis$auto$design, "panel_missing"),
          identical(pm_dis$periods_used, c("2", "3")),
          isTRUE(same_fit(pm_dir, pm_dis, skip = c("call", "auto", "note"))),
          all(pm_dir$note %in% pm_dis$note),
          any(grepl("was given a `time =` index", pm_dis$note, fixed = TRUE)))
## the diagnosis says so too, and one within-cell index is required, not two
chk_st <- mwperm_check(index = c("i", "j"), time = "t", data = st,
                       design = "irregular")
stopifnot(identical(chk_st$design, "panel_missing"),
          identical(chk_st$roles$time, "t"))
expect_err(mwperm_check(index = c("i", "j"), time = "t", rep = "t",
                        data = st, design = "irregular"),
           "not both")

passed("test-main.R")
