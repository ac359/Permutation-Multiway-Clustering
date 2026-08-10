## Detection-safety tests: assumption forks and the time-role assignment must
## never pick a potentially anti-conservative branch silently.
## Base-R stopifnot style; fast.
library(mwperm)

warns_of <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  w
}

## shared fixture: complete balanced 6 x 6 x 4 crossed array
set.seed(42)
g3 <- expand.grid(a = 1:6, b = 1:6, tt = 1:4)
N  <- nrow(g3)
yv <- rnorm(N)
dv <- rnorm(N)

## ---- 1. F5.5: a name-alone time assignment must warn -----------------------
## 'period' is a genuine cluster (6 levels); the true time dimension is 'g'
## (4 levels, strictly fewest -- strongly time-like by value). The name match
## must still win the role (dispatch behaviour is frozen) but now carries a
## warning with override instructions.
idx_f55 <- list(period = g3$a, i = g3$b, g = g3$tt)
chk <- mwperm_check(index = idx_f55)
stopifnot(identical(chk$design, "panel"),
          identical(chk$roles$time, "period"))     # assignment unchanged
stopifnot(any(grepl("name", chk$warnings, ignore.case = TRUE) &
              grepl("time =", chk$warnings, fixed = TRUE)))

w1 <- warns_of(fit1 <- mwperm(y = yv, d = dv, index = idx_f55,
                              conf_int = FALSE, seed = 1, verbose = FALSE))
stopifnot(inherits(fit1, "mwperm"),
          any(grepl("name", w1, ignore.case = TRUE)),
          any(grepl("name", fit1$note,
                    ignore.case = TRUE)))  # kept on the object

## ---- 2. a values-corroborated name assignment stays silent -----------------
## 'year' is regularly spaced with strictly fewer levels than both clusters:
## the values single it out, so no warning fires (the canonical panel path
## must not become noisy).
idx_ok <- list(i = g3$a, j = g3$b, year = 2000L + g3$tt)
chk2 <- mwperm_check(index = idx_ok)
stopifnot(identical(chk2$design, "panel"),
          identical(chk2$roles$time, "year"),
          length(chk2$warnings) == 0L)
w2 <- warns_of(mwperm(y = yv, d = dv, index = idx_ok,
                      conf_int = FALSE, seed = 1, verbose = FALSE))
stopifnot(length(w2) == 0L)

## ---- 3. F3.4: forcing threeway over a time-like index must warn ------------
chk3 <- mwperm_check(index = idx_ok, design = "threeway")
stopifnot(identical(chk3$design, "threeway"),
          any(grepl("time-like", chk3$warnings) &
              grepl("panel", chk3$warnings)))
w3 <- warns_of(fit3 <- mwperm(y = yv, d = dv, index = idx_ok,
                              design = "threeway",
                              conf_int = FALSE, seed = 1, verbose = FALSE))
stopifnot(inherits(fit3, "mwperm"), any(grepl("time-like", w3)))

## forced threeway with no time-like evidence stays silent
g5 <- expand.grid(a = 1:5, b = 1:5,
                  cc = 1:5)   # equal levels: nothing stands out
set.seed(7)
y5 <- rnorm(nrow(g5))
d5 <- rnorm(nrow(g5))
idx_ex <- list(a = g5$a, b = g5$b, cc = g5$cc)
chk4 <- mwperm_check(index = idx_ex, design = "threeway")
stopifnot(length(chk4$warnings) == 0L)
w4 <- warns_of(mwperm(y = y5, d = d5, index = idx_ex, design = "threeway",
                      conf_int = FALSE, seed = 1, verbose = FALSE))
stopifnot(length(w4) == 0L)

## ---- 4. ambiguous fork: default is announced and worded correctly ----------
## No name evidence and non-discriminating values -> panel default with the
## third index held fixed; the warning must say the default protects only if
## the PERMUTED pair is exchangeable (not "valid even if all three are").
idx_amb <- list(a = g3$a, b = g3$b, cc = g3$tt)
chk5 <- mwperm_check(index = idx_amb)
stopifnot(identical(chk5$design, "panel"),
          identical(chk5$roles$time, "cc"),
          any(grepl("permuted", chk5$warnings) &
              grepl("time =", chk5$warnings, fixed = TRUE)))

## ---- 5. layout fork warning is raised at fit time --------------------------
g2 <- expand.grid(i = 1:5, j = 1:5)
g2 <- rbind(g2, g2, g2)                          # three replicates per cell
set.seed(11)
yl <- rnorm(nrow(g2))
dl <- rnorm(nrow(g2))
w5 <- warns_of(fit5 <- mwperm(y = yl, d = dl,
                              index = list(i = g2$i, j = g2$j),
                              conf_int = FALSE, seed = 1, verbose = FALSE))
stopifnot(inherits(fit5, "mwperm"),
          any(grepl("replication", w5)))

cat("test-detection.R: all assertions passed\n")
