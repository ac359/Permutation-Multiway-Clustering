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
expect_err(mwperm_check(index = list(i = gm$i, j = gm$j),
                        time = rep(1, nrow(gm))), "complete")

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

passed("test-main.R")
