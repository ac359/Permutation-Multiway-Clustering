## Tests for the plot.mwperm figure types and mwperm_save. Base-R stopifnot
## style; everything draws on a null device so nothing hits disk except the
## mwperm_save files (written to tempdir). Fast: small designs, few reps.
library(mwperm)

msg_of <- function(expr) tryCatch({expr; NA_character_},
                                  error = function(e) conditionMessage(e))
warn_of <- function(expr) {                 # collect warning() output, run expr
  out <- character(0)
  withCallingHandlers(expr, warning = function(w) {
    out <<- c(out, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  out
}
msgs_of <- function(expr) {                 # collect message() output, run expr
  out <- character(0)
  withCallingHandlers(expr, message = function(m) {
    out <<- c(out, conditionMessage(m)); invokeRestart("muffleMessage")
  })
  out
}
## reach a package internal both when installed (namespace) and when the R/
## files are source()d directly for development (global env)
internal <- function(nm) {
  if (exists(nm, mode = "function", inherits = TRUE))
    get(nm, mode = "function", inherits = TRUE)
  else get(nm, envir = asNamespace("mwperm"))
}

grDevices::pdf(NULL)                        # null device: no files written

## ---- shared toy data --------------------------------------------------------
set.seed(42)
n <- 21L                                    # >= 20 clusters/dim so a 95% CI exists
g <- expand.grid(i = seq_len(n), j = seq_len(n)); N <- nrow(g)
d1 <- rnorm(n)[g$i] + rnorm(N)
y1 <- rnorm(n)[g$i] + rnorm(n)[g$j] + 0.4 * d1 + rnorm(N)

## ---- 1. d = 1 with CI: every type runs, invisible return, par restored -----
fit1 <- mwperm_dyadic(y1, d1, row = g$i, col = g$j, seed = 3, n_reps = 7)
stopifnot(!is.null(fit1$conf_int))
vis <- withVisible(plot(fit1))
stopifnot(!vis$visible, identical(vis$value, fit1))
mar0 <- graphics::par("mar"); mfrow0 <- graphics::par("mfrow")
for (tp in c("auto", "coef", "stability", "all")) plot(fit1, type = tp)
stopifnot(identical(graphics::par("mar"), mar0),
          identical(graphics::par("mfrow"), mfrow0))   # S6: par() restored

## single rep (the default call): coef figure by default, strip diagnostic
fit1b <- mwperm_dyadic(y1, d1, row = g$i, col = g$j, seed = 9)
plot(fit1b)
plot(fit1b, type = "stability")

## ---- 2. d = 2: forest (coef), joint region, all ------------------------------
D2 <- cbind(a = d1, b = rnorm(n)[g$j] + rnorm(N))
y2 <- y1 - 0.2 * D2[, "b"]
fit2 <- mwperm_dyadic(y2, D2, row = g$i, col = g$j, seed = 4,
                      grid = list(seq(-0.5, 1.2, length.out = 11),
                                  seq(-1.0, 0.6, length.out = 11)))
stopifnot(!is.null(fit2$conf_region), nrow(fit2$conf_region) >= 1L)
for (tp in c("auto", "coef", "region", "stability", "all")) plot(fit2, type = tp)

## vector null: the forest draws per-row null crosses instead of one line
fit2b <- mwperm_dyadic(y2, D2, row = g$i, col = g$j, seed = 4,
                       beta_null = c(0.4, -0.2),
                       grid = list(seq(-0.5, 1.2, length.out = 11),
                                   seq(-1.0, 0.6, length.out = 11)))
plot(fit2b, type = "coef")

## empty region (grid far from the acceptance set): fall back, never error
fit2e <- mwperm_dyadic(y2, D2, row = g$i, col = g$j, seed = 8,
                       grid = list(c(100, 101), c(100, 101)))
plot(fit2e)
stopifnot(length(msgs_of(plot(fit2e, type = "region"))) >= 1L)

## ---- 3. d = 3: forest from conf_box; "region" falls back (needs d = 2) ------
D3 <- cbind(D2, c = rnorm(N))
y3 <- y2 + 0.1 * D3[, "c"]
fit3 <- mwperm_dyadic(y3, D3, row = g$i, col = g$j, seed = 5,
                      grid = list(seq(0.0, 0.8, length.out = 5),
                                  seq(-0.6, 0.2, length.out = 5),
                                  seq(-0.3, 0.5, length.out = 5)))
plot(fit3)
stopifnot(any(grepl("two coefficients", msgs_of(plot(fit3, type = "region")))))

## ---- 4. no confidence set stored: message-and-fallback ----------------------
fit_no <- mwperm_dyadic(y1, d1, row = g$i, col = g$j, seed = 6, conf_int = FALSE)
stopifnot(any(grepl("stability", msgs_of(plot(fit_no, type = "coef")))))
plot(fit_no)                                # auto -> stability, silently fine

## coarse resolution (6 x 6 -> 1/(K+1) > alpha, conf_int NULL with a note)
g6 <- expand.grid(i = 1:6, j = 1:6)
fit6 <- mwperm_dyadic(rnorm(36), rnorm(36), row = g6$i, col = g6$j, seed = 7)
for (tp in c("auto", "coef", "stability", "all")) plot(fit6, type = tp)

## ---- 5. reserved types message and fall back, never error -------------------
stopifnot(length(msgs_of(plot(fit1, type = "null"))) >= 1L)
stopifnot(length(msgs_of(plot(fit1, type = "profile"))) >= 1L)

## ---- 6. every other design draws without error ------------------------------
gp <- expand.grid(i = 1:6, j = 1:6, t = 1:3); Np <- nrow(gp)
dp <- rnorm(6)[gp$i] + rnorm(Np)
yp <- rnorm(6)[gp$i] + rnorm(6)[gp$j] + c(0, 1, 2)[gp$t] + 0.5 * dp + rnorm(Np)
fitp <- mwperm_panel(yp, dp, row = gp$i, col = gp$j, time = gp$t,
                     seed = 12, conf_int = FALSE)

gt <- expand.grid(a = 1:5, b = 1:5, c = 1:4); Nt <- nrow(gt)
dt_ <- rnorm(Nt)
yt <- rnorm(5)[gt$a] + rnorm(5)[gt$b] + rnorm(4)[gt$c] + 0.3 * dt_ + rnorm(Nt)
fitt <- mwperm_threeway(yt, dt_, id1 = gt$a, id2 = gt$b, id3 = gt$c,
                        seed = 13, conf_int = FALSE)

gl <- expand.grid(i = 1:5, j = 1:5, l = 1:4); Nl <- nrow(gl)
dl <- rnorm(Nl)
eta <- rnorm(25L)[(gl$i - 1L) * 5L + gl$j]
yl <- eta + 0.3 * dl + rnorm(Nl)
fitl <- mwperm_layout(yl, dl, row = gl$i, col = gl$j, seed = 14,
                      conf_int = FALSE)

gm <- expand.grid(i = 1:10, j = 1:10); gm <- gm[gm$i != gm$j, ]  # drop diagonal
Nm <- nrow(gm)
dm <- rnorm(10)[gm$i] + rnorm(Nm)
ym <- rnorm(10)[gm$i] + rnorm(10)[gm$j] + 0.5 * dm + rnorm(Nm)
fitm <- mwperm_missing(ym, dm, row = gm$i, col = gm$j, min_block = 3,
                       seed = 11, conf_int = FALSE)

for (f in list(fitp, fitt, fitl, fitm)) {
  out <- withVisible(plot(f))
  stopifnot(!out$visible, identical(out$value, f))
  plot(f, type = "all")
}

## ---- 7. style layer: Okabe-Ito defaults, S1 (never colour alone), overrides -
sty <- internal(".mwperm_style")
s <- sty()
stopifnot(identical(s$col_estimate, "#0072B2"),      # Okabe-Ito blue
          identical(s$col_null, "#D55E00"))          # Okabe-Ito vermillion
stopifnot(s$pch_estimate != s$pch_null,              # marks differ, not just colour
          s$lty_null != s$lty_alpha)                 # reference lines differ by lty
s2 <- sty(col_estimate = "black", lwd_interval = 3)
stopifnot(identical(s2$col_estimate, "black"), s2$lwd_interval == 3)
m <- msg_of(sty(no_such_element = 1))
stopifnot(grepl("no_such_element", m, fixed = TRUE))

## overrides through plot(...) run; explicit graphics args no longer crash
plot(fit1, col_estimate = "black", lwd_interval = 3)
plot(fit1, main = "custom title", sub = "", xlab = "beta")
w <- warn_of(plot(fit1, bogus_argument = 1))         # unknown ...: warn + draw
stopifnot(any(grepl("bogus_argument", w, fixed = TRUE)))

## ---- 8. mwperm_save: journal-dimension export -------------------------------
tf <- file.path(tempdir(), "mwperm-test-fig.png")
out <- mwperm_save(fit1, tf, width = "single")
stopifnot(identical(out, tf), file.exists(tf), file.size(tf) > 0)
unlink(tf)
tfp <- file.path(tempdir(), "mwperm-test-fig.pdf")
mwperm_save(fit1, tfp, width = "double", type = "stability")
stopifnot(file.exists(tfp)); unlink(tfp)
m <- msg_of(mwperm_save(fit1, file.path(tempdir(), "fig.bmp")))
stopifnot(grepl("extension", m, fixed = TRUE))
m <- msg_of(mwperm_save(list(), file.path(tempdir(), "fig.png")))
stopifnot(grepl("mwperm", m, fixed = TRUE))
m <- msg_of(mwperm_save(fit1, file.path(tempdir(), "fig.png"), width = -1))
stopifnot(grepl("`width`", m, fixed = TRUE))

grDevices::dev.off()
cat("test-plot.R: all assertions passed\n")
