## Code-path regression tests for branches the design-focused files miss: the
## design-diagnosis printer, forced-design validation, name resolution, the
## L0 balancing path, explicit-grid confidence sets, the constructed
## disconnected-acceptance island guard, and the figure-export devices.
## Complements the behaviour suites; added to bring shipped-test coverage of
## R/ above 90%. Base-R stopifnot style; fast.
library(mwperm)

msg_of <- function(expr)
  tryCatch({
    expr
    NA_character_
  }, error = function(e) conditionMessage(e))
expect_err <- function(expr, pattern) {
  m <- msg_of(expr)
  stopifnot(!is.na(m), grepl(pattern, m, fixed = TRUE))
}
warns_of <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  w
}

## ---- fixtures ---------------------------------------------------------------
set.seed(202)
g2 <- expand.grid(i = 1:6, j = 1:6)
y2 <- rnorm(36)
d2 <- rnorm(36)
df2 <- data.frame(i = g2$i, j = g2$j, yy = y2, dd = d2, xx = rnorm(36))
g3 <- expand.grid(i = 1:5, j = 1:5, year = 1:3)
y3 <- rnorm(nrow(g3))
d3 <- rnorm(nrow(g3))
gm <- g2[g2$i != g2$j, ]
ym <- rnorm(nrow(gm))
dm <- rnorm(nrow(gm))

## ---- 1. the design-diagnosis printer ----------------------------------------
out_d <- paste(capture.output(print(
  mwperm_check(index = list(i = g2$i, j = g2$j)))), collapse = "\n")
stopifnot(grepl("dyadic", out_d), grepl("Resolution", out_d),
          grepl("TOO COARSE", out_d))            # 6 clusters < 20
out_m <- paste(capture.output(print(
  mwperm_check(index = list(i = gm$i, j = gm$j)))), collapse = "\n")
stopifnot(grepl("missing", out_m), grepl("biclique", out_m))
out_p <- paste(capture.output(print(
  mwperm_check(index = list(i = g3$i, j = g3$j, year = g3$year)))),
  collapse = "\n")
stopifnot(grepl("panel", out_p))
## a diagnosis with warnings renders them with the "!" prefix
out_w <- paste(capture.output(print(
  mwperm_check(index = list(period = g3$i, a = g3$j, b = g3$year)))),
  collapse = "\n")
stopifnot(grepl("!", out_w, fixed = TRUE))

## ---- 2. index resolution + forced-design validation -------------------------
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
gm3 <- g3[-1, ]                                   # incomplete 3-index array
expect_err(mwperm_check(index = list(i = gm3$i, j = gm3$j, t = gm3$year),
                        design = "threeway"), "complete balanced")
expect_err(mwperm_check(index = list(i = g3$i, j = g3$j, t = g3$year),
                        design = "layout"), "exactly 2")
expect_err(mwperm_check(index = list(i = g3$i, j = g3$j, t = g3$year),
                        time = g3$year), "Too many clustering dimensions")
expect_err(mwperm_check(index = list(i = g3$i, j = g3$j, t = g3$year),
                        rep = g3$year), "Too many clustering dimensions")
expect_err(mwperm_check(index = list(i = gm$i, j = gm$j), time = rep(1,
                                                                     nrow(gm))),
           "complete")
## forced missing on 2 complete/incomplete indices: K depends on the blocks
chk_fm <- mwperm_check(index = list(i = gm$i, j = gm$j), design = "missing")
stopifnot(identical(chk_fm$design, "missing"), is.na(chk_fm$K_default))

## mwperm(): data-column resolution incl. multi-column d, and check_arg
f_nm <- mwperm(y = "yy", d = c("dd", "xx"), index = c("i", "j"), data = df2,
               seed = 2, n_reps = 1, conf_int = FALSE, verbose = FALSE)
stopifnot(length(f_nm$estimate) == 2L,
          identical(f_nm$d_names, c("dd", "xx")))
expect_err(mwperm(y = "zzz", d = "dd", index = c("i", "j"), data = df2,
                  verbose = FALSE), "`y`")
w_arg <- warns_of(mwperm(y = "yy", d = "dd", index = c("i", "j"), data = df2,
                         time_fe = FALSE, L0 = 3, seed = 2, n_reps = 1,
                         conf_int = FALSE, verbose = FALSE))
stopifnot(any(grepl("`time_fe`", w_arg, fixed = TRUE)),
          any(grepl("`L0`", w_arg, fixed = TRUE)))

## ---- 3. layout: L0 balancing and the no-power warning -----------------------
gl <- expand.grid(l = 1:6, i = 1:5, j = 1:5)
gl <- gl[!(gl$l > 3 & gl$i == 1), ]               # unequal cell sizes
yl <- rnorm(nrow(gl))
dl <- rnorm(nrow(gl))
fL  <- mwperm_layout(yl, dl, row = gl$i, col = gl$j, L0 = 3, seed = 2,
                     n_reps = 2, conf_int = FALSE)
fL2 <- mwperm_layout(yl, dl, row = gl$i, col = gl$j, L0 = 3, seed = 2,
                     n_reps = 2, conf_int = FALSE)
stopifnot(any(grepl("Balanced to L0", fL$note)),
          fL$K <= 2L,                             # K capped at L0 - 1
          identical(fL$pvalue, fL2$pvalue))       # L0 downsampling is seeded
dcc <- rnorm(25)[(gl$i - 1L) * 5L + gl$j]         # constant within every cell
w_cc <- warns_of(mwperm_layout(yl, dcc, row = gl$i, col = gl$j, seed = 1,
                               n_reps = 1, conf_int = FALSE))
stopifnot(any(grepl("constant within every cell", w_cc)))

## ---- 4. explicit-grid and joint-region confidence paths ---------------------
## (alpha = 0.4 keeps 1/(K+1) = 1/6 attainable on these 6x6 fixtures)
## NB grid mode now inverts the MEDIAN p-value across reps (not the old union /
## maximum); the median-vs-union contrast and the edge/NA guards are covered in
## detail by test-invert-ci-grid.R (NEWS 0.2.0 "Bug fix (changes intervals)").
f_g1 <- mwperm_dyadic(y2, d2, row = g2$i, col = g2$j, seed = 1, n_reps = 2,
                      alpha = 0.4, grid = seq(-3, 3, by = 0.25))
stopifnot(!is.null(f_g1$conf_int), all(is.finite(f_g1$conf_int)))
## d = 4: the 21^4 default region grid is refused with a note, never computed
f_d4 <- mwperm_dyadic(y2, matrix(rnorm(36 * 4), 36), row = g2$i, col = g2$j,
                      seed = 1, n_reps = 1, alpha = 0.4)
stopifnot(any(grepl("too large", f_d4$note)))
## d = 2 with an explicit per-coefficient grid: region + box computed
f_d2 <- mwperm_dyadic(y2, cbind(a = d2, b = rnorm(36)), row = g2$i,
                      col = g2$j, seed = 1, n_reps = 1, alpha = 0.4,
                      grid = list(seq(-4, 4, by = 0.5), seq(-4, 4, by = 0.5)))
stopifnot(!is.null(f_d2$conf_box) || any(grepl("region", f_d2$note)))
## d = 2 with too many explicit grid points: refused with a note
f_d2b <- mwperm_dyadic(y2, cbind(a = d2, b = rnorm(36)), row = g2$i,
                       col = g2$j, seed = 1, n_reps = 1, alpha = 0.4,
                       grid = list(seq_len(150), seq_len(150)))
stopifnot(any(grepl("points", f_d2b$note)))
## degenerate d with conf_int = TRUE: CI search falls back and stays honest
w_dg <- warns_of(f_dg <- mwperm_dyadic(y2, rep(1, 36), row = g2$i, col = g2$j,
                                       seed = 1, n_reps = 1, alpha = 0.4))
stopifnot(any(grepl("`d`", w_dg, fixed = TRUE)),
          f_dg$pvalue == 1, all(is.infinite(f_dg$conf_int)))

## ---- 5. island guard: constructed disconnected acceptance set ---------------
## a_k(b) = |t_k - b|,
## b_k(b) = |b|/2 with t = (-5 x10, +5 x9) gives the exact acceptance set
## [-10, -10/3] U [10/3, 10] at alpha = .05; the reported interval must be
## the HULL spanning both islands, flagged as disconnected.
Ki <- 19L
prep_i <- list(u = matrix(c(rep(-5, 10), rep(5, 9)), 1),
               v = matrix(0, 1, Ki),
               M = array(1, c(1, 1, Ki)), W = array(0.5, c(1, 1, Ki)),
               K = Ki, Kp1 = 20L, d = 1L, has_perm_D = TRUE)
ci_i <- mwperm:::.invert_ci(list(prep_i), alpha = 0.05, centre = 0, scale = 1,
                            y = c(0, 1), D = matrix(c(0, 1)))
stopifnot(isTRUE(attr(ci_i, "disconnected")),
          ci_i[1] <= -10 + 0.01, ci_i[2] >= 10 - 0.01)

## ---- 6. internal one-liners -------------------------------------------------
stopifnot(mwperm:::.fmt_p(NA) == "NA",
          mwperm:::.fmt_p(1e-5) == "< 0.001")
expect_err(mwperm:::.build_obs_perms(cbind(1:4, 1:4), list(NULL, NULL)),
           "At least one dimension")
expect_err(mwperm:::.plapply(as.list(1:3), identity, n_cores = NA),
           "`n_cores`")
gp <- expand.grid(i = 1:4, j = 1:4, t = 1:2)
expect_err(mwperm_panel(c(rnorm(nrow(gp)), 1), rnorm(nrow(gp) + 1L),
                        row = c(gp$i, 1), col = c(gp$j, 1),
                        time = c(gp$t, 1)), "at most once")

## ---- 7. figure export: every device + argument validation -------------------
fitp <- mwperm_dyadic(y2, d2, row = g2$i, col = g2$j, seed = 1, n_reps = 6,
                      alpha = 0.4)
for (ext in c("pdf", "png", "tiff", "jpeg")) {
  f <- file.path(tempdir(), paste0("mwperm-cov.", ext))
  mwperm_save(fitp, f)
  stopifnot(file.exists(f), file.size(f) > 0)
  unlink(f)
}
mwperm_save(fitp, file.path(tempdir(), "mwperm-all.pdf"), width = "double",
            type = "all", pointsize = 9)
unlink(file.path(tempdir(), "mwperm-all.pdf"))
expect_err(mwperm_save(fitp, 1), "`file`")
expect_err(mwperm_save(fitp, "a.png", width = -1), "`width`")
expect_err(mwperm_save(fitp, "a.png", height = 0), "`height`")
expect_err(mwperm_save(fitp, "a.png", res = 10), "`res`")
expect_err(mwperm_save("x", "a.png"), "mwperm")
## the dot-strip stability variant (n_reps < 5) draws without error
f2r <- mwperm_dyadic(y2, d2, row = g2$i, col = g2$j, seed = 1, n_reps = 2,
                     conf_int = FALSE)
grDevices::pdf(NULL)
plot(f2r, type = "stability")
grDevices::dev.off()

cat("test-paths.R: all assertions passed\n")
