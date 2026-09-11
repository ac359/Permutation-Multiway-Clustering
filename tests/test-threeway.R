## mwperm_threeway() -- the three-way design, condition InvA in three
## dimensions.
##
## All three cluster dimensions are permuted jointly by a common group order
## (Guo, Toulis & Wang 2026, Section 6.1), which is exact under separate
## exchangeability in all three. The array must be complete and balanced: a
## joint relabelling of three dimensions is only defined when every cell of the
## id1 x id2 x id3 grid is observed exactly once.
##
## Cross-cutting contracts live elsewhere; see tests/README.md. The warning
## raised when a time-like index is forced into this design is in test-main.R,
## because the choice is made by the detector, not here.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))

## ---- fixture: a 5 x 6 x 7 array with three-way random effects -------------
set.seed(2)
g <- expand.grid(a = 1:5, b = 1:6, c = 1:7)
N <- nrow(g)
d <- rnorm(N)
eps <- rnorm(5)[g$a] + rnorm(6)[g$b] + rnorm(7)[g$c] + rnorm(N)
y <- 0.8 * d + eps

## ---- 1. structure, seeded reproducibility, and K --------------------------
## Default K = min(n1, n2, n3) - 1: the joint group order cannot exceed the
## smallest dimension, because every dimension is relabelled by it.
fit <- mwperm_threeway(y, d, id1 = g$a, id2 = g$b, id3 = g$c, seed = 1,
                       conf_int = FALSE)
stopifnot(inherits(fit, "mwperm"),
          identical(fit$type, "threeway"),
          fit$K == 4L, fit$n_perm == 5L,
          fit$n_obs == N,
          identical(unname(fit$n_clusters), c(5L, 6L, 7L)),
          identical(names(fit$n_clusters), c("id1", "id2", "id3")))
fit_b <- mwperm_threeway(y, d, id1 = g$a, id2 = g$b, id3 = g$c, seed = 1,
                         conf_int = FALSE)
stopifnot(isTRUE(same_fit(fit, fit_b)))
expect_err(mwperm_threeway(y, d, id1 = g$a, id2 = g$b, id3 = g$c, K = 5L),
           "K + 1")

## ---- 2. the p-value responds to signal ------------------------------------
## A strong true effect drives the p-value to its floor 1/(K+1); the same
## design with the outcome carrying no d at all does not reject. Size and
## power are measured properly by the Monte-Carlo suite in inst/replication/ --
## this is the sanity check that the front end is wired to the right column.
strong <- mwperm_threeway(eps + 5 * d, d, id1 = g$a, id2 = g$b, id3 = g$c,
                          seed = 1, conf_int = FALSE)
stopifnot(identical(strong$pvalue, strong$p_floor))
none <- mwperm_threeway(eps, d, id1 = g$a, id2 = g$b, id3 = g$c, seed = 1,
                        conf_int = FALSE)
stopifnot(none$pvalue > 0.05)

## ---- 3. row-order invariance ----------------------------------------------
o <- sample(N)
fit_o <- mwperm_threeway(y[o], d[o], id1 = g$a[o], id2 = g$b[o], id3 = g$c[o],
                         seed = 1, conf_int = FALSE)
stopifnot(identical(fit$pvalue, fit_o$pvalue),
          abs(fit$estimate - fit_o$estimate) < 1e-10)

## ---- 4. the array must be complete and balanced ---------------------------
expect_err(mwperm_threeway(y[-1], d[-1], id1 = g$a[-1], id2 = g$b[-1],
                           id3 = g$c[-1]),
           "complete balanced array")
expect_err(mwperm_threeway(c(y, 1), c(d, 1), id1 = c(g$a, 1), id2 = c(g$b, 1),
                           id3 = c(g$c, 1)),
           "at most once")

## ---- 5. a joint null over several columns of d ----------------------------
## d > 1 tests H0: beta = b jointly; the acceptance region is a set of vectors,
## so conf_int is NULL and the region/box carry the answer.
D2 <- cbind(a = d, b = rnorm(N))
fit2 <- mwperm_threeway(y + 0.3 * D2[, "b"], D2, id1 = g$a, id2 = g$b,
                        id3 = g$c, seed = 1, alpha = 0.4)
stopifnot(length(fit2$estimate) == 2L,
          identical(fit2$d_names, c("a", "b")),
          is.null(fit2$conf_int),
          !is.null(fit2$conf_box) || any(grepl("region", fit2$note)))

passed("test-threeway.R")
