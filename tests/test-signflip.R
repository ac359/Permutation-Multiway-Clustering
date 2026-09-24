## mwperm_dyadic_het() and build_flip_set() -- the sign-flip (Rademacher)
## invariant test, IPT-Het.
##
## The permutation designs need the error array to be exchangeable under
## relabelling of the clusters; that fails under heteroskedasticity. The
## sign-flip design replaces the permutation group by the group of joint
## row-and-column sign changes, valid under symmetry of the errors about zero
## and indifferent to their variances (Guo, Toulis & Wang 2026, Section 2:
## the argument holds for any invariance group). Every group element is a
## signed gather `list(g = NULL, s = +/-1)` applied through `.apply_op()`,
## which is the one generalisation the engine needed.
##
## Section 1 is the test that matters most: agreement with a corrected port of
## the author's research implementation (tests/helpers/signflip-reference.R).
## Every section states the failure mode it guards against, not the assertion.
##
##   1. corrected-reference agreement   -- the package IS Procedure 1 under
##                                         the sign-flip group
##   2. the group                       -- closure, identity once, 2^(n_flip-1)
##                                         distinct elements, all groups used
##   3. .apply_op()                     -- the widened engine contract
##   4. size under heteroskedasticity   -- validity of the new test (small MC)
##   5. resolution guard                -- the floor is 1/2^(n_flip-1), the
##                                         notes speak of n_flip, n_flip = 6
##                                         is the first with a 95% set
##   6. seed hygiene                    -- reproducible, and the caller's RNG
##                                         stream is untouched
##   7. engine properties               -- affine-in-b, parallel == serial,
##                                         every S3 method runs
##   8. the dispatcher                  -- opt-in only, dispatch identity,
##                                         argument cross-checks
##   9. permutation designs untouched   -- seeded anchor and note wording
##
## Run from the package root after installing:
##   R CMD INSTALL . && Rscript tests/test-signflip.R
## Slower, wider checks (the author's original script verbatim, from-scratch
## projectors, a cross-version bit-identity battery) are development-only
## and are not shipped with the package.
library(mwperm)
source(if (file.exists("helpers/assertions.R")) "helpers/assertions.R"
       else file.path("tests", "helpers", "assertions.R"))
source(if (file.exists("helpers/signflip-reference.R"))
         "helpers/signflip-reference.R"
       else file.path("tests", "helpers", "signflip-reference.R"))

sub_seed <- internal(".sub_seed")
apply_op <- internal(".apply_op")

## ---- 1. corrected-reference agreement -------------------------------------
## Guards against the package computing SOME sign-flip statistic rather than
## Procedure 1 with the sign-flip group: the reference shares no code with R/
## (explicit orthonormal complement via a complete QR, array layout, the
## author's own a/b formulas), differs from the author's script only in the
## two documented corrections, and must return the IDENTICAL p-value when
## handed the same flip-group assignment. A last-bit disagreement in a_k vs
## b_k would show up here as a moved p-value.
set.seed(2026)
n12 <- 12L
nf5 <- 5L
dat <- sim_gravity_het(n12, b = 0, rho1 = 0.5, rho2 = 0.3)
Y_arr <- array(dat$y, dim = c(n12, n12, 1L))
D_arr <- array(dat$d, dim = c(n12, n12, 1L))
X_arr <- array(cbind(1, dat$x), dim = c(n12, n12, 1L, 3L))
for (sd in 1:8) for (b0 in c(0, 0.25)) {
  fit <- mwperm_dyadic_het(dat$y, dat$d, x = dat$x, row = dat$i, col = dat$j,
                           n_flip = nf5, n_reps = 1, seed = sd,
                           beta_null = b0, conf_int = FALSE)
  ## the assignment the fit drew: rep 1 of seed `sd` is sub-seed (sd, 1)
  F <- build_flip_set(n12, n12, nf5, seed = sub_seed(sd, 1L))
  p_ref <- rpt_signflip_ref(Y_arr, X_arr, D_arr,
                            attr(F, "row_groups"), attr(F, "col_groups"),
                            nf5, beta_null = b0)
  stopifnot(identical(fit$pvalue, p_ref))
}

## ---- 2. the group -----------------------------------------------------------
## Guards against a set that is not a group (Theorem 1 needs closure), the
## identity appearing twice (the {s, -s} kernel: 2^n_flip enumeration would
## list every element twice), or a wrong count (the resolution 1/2^(n_flip-1)
## is only right if there are exactly that many DISTINCT elements).
F <- build_flip_set(n_row = 7, n_col = 5, n_flip = 4, seed = 3)
stopifnot(length(F) == 2^3,
          identical(attr(F, "group_order"), 8L),
          identical(attr(F, "n_flip"), 4L),
          length(attr(F, "row_groups")) == 7L,
          length(attr(F, "col_groups")) == 5L,
          ## every flip group is used (the unused-group guard)
          length(unique(c(attr(F, "row_groups"), attr(F, "col_groups")))) == 4L,
          all(F[[1]]$row == 1), all(F[[1]]$col == 1),        # identity first
          all(vapply(F, function(e) all(abs(e$row) == 1) && all(abs(e$col) == 1),
                     logical(1))))
## observation-level action of every element, as an m x n sign matrix
act <- lapply(F, function(e) outer(e$row, e$col))
key <- vapply(act, function(a) paste(a, collapse = ""), character(1))
stopifnot(!anyDuplicated(key),                              # all distinct
          sum(vapply(act, function(a) all(a == 1), logical(1))) == 1L)
## closure: the product of any two elements is an element
for (a in seq_along(F)) for (b in seq_along(F)) {
  prod <- paste(act[[a]] * act[[b]], collapse = "")
  stopifnot(prod %in% key)
}
## inverses: every element is its own inverse (a sign flip applied twice)
stopifnot(all(vapply(act, function(a) all(a * a == 1), logical(1))))
## the coset representatives pin coordinate 1 to +1, so no s and -s both appear
S <- attr(F, "signs")
stopifnot(all(S[, 1] == 1), nrow(S) == 8L, ncol(S) == 4L,
          !anyDuplicated(S))

## argument errors name the argument
expect_err(build_flip_set(4, 4, 1), "`n_flip`")
expect_err(build_flip_set(3, 2, 6), "`n_flip` = 6 exceeds")
expect_err(build_flip_set(0, 4, 2), "`n_row`")

## ---- 3. .apply_op(): the signed-gather contract -----------------------------
## Guards the one engine change: a bare integer vector must be applied by
## EXACTLY the indexing the permutation front ends always used (bit-identical
## seeded output), and the two slots must mean gather-then-flip with NULL
## as the identity in each.
M <- matrix(rnorm(12), 4, 3)
v <- rnorm(4)
g <- c(3L, 1L, 4L, 2L)
s <- c(1, -1, -1, 1)
stopifnot(identical(apply_op(g, M), M[g, , drop = FALSE]),
          identical(apply_op(g, v), v[g]),
          identical(apply_op(list(g = g, s = NULL), M), M[g, , drop = FALSE]),
          identical(apply_op(list(g = NULL, s = s), M), M * s),
          identical(apply_op(list(g = NULL, s = s), v), v * s),
          identical(apply_op(list(g = g, s = s), M), M[g, , drop = FALSE] * s),
          identical(apply_op(list(g = NULL, s = NULL), M), M))

## ---- 4. size under heteroskedasticity ---------------------------------------
## Guards the reason the test exists: on the heteroskedastic gravity DGP (the
## author's simulate_gravity_model(), error sd increasing in the gravity mean
## and in the distance covariate) the sign-flip test must reject a true null
## at or below the nominal rate, heteroskedastic or not. Small and seeded (R
## CMD check runs this), so it asserts validity of the NEW test only, with
## Monte-Carlo slack; the permutation test's over-rejection on the same DGP
## -- 0.06-0.09 against the sign-flip test's 0.02 at n = 25 and 300
## replications, rising with the strength of the heteroskedasticity -- is a
## comparison between two different p-value grids that 40 draws at a small n
## cannot resolve, and is recorded in the design note instead.
n_mc <- 20L
R_mc <- 40L
alpha_mc <- 0.10                      # order 32: p <= 0.10 is attainable
size_flip <- function(rho1, rho2) {
  p <- numeric(R_mc)
  for (r in seq_len(R_mc)) {
    set.seed(500 + r)
    dd <- sim_gravity_het(n_mc, b = 0, rho1 = rho1, rho2 = rho2)
    p[r] <- mwperm_dyadic_het(dd$y, dd$d, x = dd$x, row = dd$i, col = dd$j,
                              n_flip = 6, n_reps = 1, seed = r,
                              conf_int = FALSE, alpha = alpha_mc)$pvalue
  }
  mean(p <= alpha_mc)
}
slack <- 2 * sqrt(alpha_mc * (1 - alpha_mc) / R_mc)   # two MC standard errors
stopifnot(size_flip(0.9, 0.6) <= alpha_mc + slack,   # strong heteroskedasticity
          size_flip(0.5, 0.3) <= alpha_mc + slack,   # mild
          size_flip(0.0, 0.0) <= alpha_mc + slack)   # homoskedastic

## ---- 5. resolution guard -----------------------------------------------------
## Guards the engine's gate seeing the RIGHT group order: 2^(n_flip - 1), not
## n_flip + 1. At n_flip = 5 the floor is 1/16 = 0.0625 > 0.05, so no 95% set
## exists and the note must say what the design needs in the sign-flip
## vocabulary; at n_flip = 6 the floor is 1/32 and the set is attainable.
data(trade_dyadic)
f5 <- with(trade_dyadic,
           mwperm_dyadic_het(y = log_trade, d = log_dist,
                             x = cbind(log_gdp_i, log_gdp_j),
                             row = importer, col = exporter,
                             n_flip = 5, n_reps = 2, seed = 1))
invisible(capture.output(s5 <- summary(f5)))    # summary() prints, then returns
stopifnot(is.null(f5$conf_int), identical(f5$K, 15L), identical(f5$n_perm, 16L),
          identical(f5$n_flip, 5L), identical(f5$p_floor, 1 / 16),
          identical(f5$resolution, 1 / 16),
          any(grepl("1/2^(n_flip-1) = 0.0625", f5$note, fixed = TRUE)),
          any(grepl("n_flip >= 6", f5$note, fixed = TRUE)),
          !any(grepl("smallest permuted dimension", f5$note, fixed = TRUE)),
          is.na(s5$ipt_ci_low))
expect_err(confint(f5), "n_flip >= 6")
f6 <- with(trade_dyadic,
           mwperm_dyadic_het(y = log_trade, d = log_dist,
                             x = cbind(log_gdp_i, log_gdp_j),
                             row = importer, col = exporter,
                             n_flip = 6, n_reps = 2, seed = 1))
stopifnot(length(f6$conf_int) == 2L, all(is.finite(f6$conf_int)),
          f6$conf_int[1] < f6$conf_int[2],
          identical(f6$ci_method, "exact"), identical(f6$n_perm, 32L),
          f6$conf_int[1] <= -1, f6$conf_int[2] >= -1,   # covers the truth
          identical(f6$type, "dyadic (sign-flip / heteroskedasticity-robust)"))
## The default n_flip follows the resolution rule (0.4.2): the smallest
## n_flip >= 2 whose REPORTED p-value floor -- 1/2^(n_flip - 1), doubled under
## aggregate = "median2" -- is at most alpha, i.e. 2^(n_flip - 1) >= m/alpha.
## At alpha = 0.05 that is 6 (order 32, floor 1/32) under "median" and 7
## (order 64) under "median2"; at alpha = 0.01, 8 and 9. Then capped at
## min(n_row, n_col) exactly as before.
dnf <- internal(".default_n_flip")
stopifnot(identical(dnf(NULL, 40L, 40L, alpha = 0.05, aggregate = "median"), 6L),
          identical(dnf(NULL, 40L, 40L, alpha = 0.05, aggregate = "median2"), 7L),
          identical(dnf(NULL, 40L, 40L, alpha = 0.01, aggregate = "median"), 8L),
          identical(dnf(NULL, 40L, 40L, alpha = 0.01, aggregate = "median2"), 9L),
          identical(dnf(NULL, 40L, 40L, alpha = 0.5, aggregate = "median"), 2L),
          identical(dnf(NULL, 40L, 40L, alpha = 0.2, aggregate = "median"), 4L),
          identical(dnf(NULL, 5L, 9L, alpha = 0.05, aggregate = "median"), 5L),
          ## an explicit value is honoured whatever alpha says
          identical(dnf(4L, 40L, 40L, alpha = 0.05, aggregate = "median"), 4L))
expect_err(dnf(NULL, 40L, 40L, alpha = 1, aggregate = "median"), "`alpha`")
f8 <- with(trade_dyadic,
           mwperm_dyadic_het(y = log_trade, d = log_dist,
                             x = cbind(log_gdp_i, log_gdp_j),
                             row = importer, col = exporter,
                             n_reps = 1, seed = 1, conf_int = FALSE))
stopifnot(identical(f8$n_flip, 6L), identical(f8$n_perm, 32L),
          identical(f8$K, 31L), f8$pvalue >= 1 / 32)
fm2 <- with(trade_dyadic,
            mwperm_dyadic_het(y = log_trade, d = log_dist,
                              x = cbind(log_gdp_i, log_gdp_j),
                              row = importer, col = exporter,
                              aggregate = "median2", n_reps = 1, seed = 1,
                              conf_int = FALSE))
f01 <- with(trade_dyadic,
            mwperm_dyadic_het(y = log_trade, d = log_dist,
                              x = cbind(log_gdp_i, log_gdp_j),
                              row = importer, col = exporter,
                              alpha = 0.01, n_reps = 1, seed = 1,
                              conf_int = FALSE))
stopifnot(identical(fm2$n_flip, 7L), identical(fm2$n_perm, 64L),
          identical(f01$n_flip, 8L), identical(f01$n_perm, 128L))
## at the default (n_flip = 6, n_reps = 10) the confidence set takes the
## exact path: 2 x 31^2 x 10 = 19,220 candidates, under the 2e5 budget
fdef <- with(trade_dyadic,
             mwperm_dyadic_het(y = log_trade, d = log_dist,
                               x = cbind(log_gdp_i, log_gdp_j),
                               row = importer, col = exporter, seed = 1))
stopifnot(identical(fdef$n_flip, 6L), identical(fdef$n_reps, 10L),
          identical(fdef$ci_method, "exact"), length(fdef$conf_int) == 2L)
## the default is capped by the smaller dimension, and the cap is enforced
g5 <- expand.grid(i = 1:5, j = 1:9)
set.seed(4)
y5 <- rnorm(45); d5 <- rnorm(45)
fc <- mwperm_dyadic_het(y5, d5, row = g5$i, col = g5$j, seed = 1, n_reps = 1,
                        conf_int = FALSE)
stopifnot(identical(fc$n_flip, 5L))
expect_err(mwperm_dyadic_het(y5, d5, row = g5$i, col = g5$j, n_flip = 6),
           "exceeds the smaller clustering dimension")
expect_err(mwperm_dyadic_het(y5, d5, row = g5$i, col = g5$j, n_flip = 1),
           "`n_flip`")
## above the cap the message names the projection count the request implies
m21 <- msg_of(with(trade_dyadic,
                   mwperm_dyadic_het(y = log_trade, d = log_dist,
                                     row = importer, col = exporter,
                                     n_flip = 21)))
stopifnot(!is.na(m21), grepl("1,048,575", m21, fixed = TRUE),
          grepl("n_flip <= 20", m21, fixed = TRUE))

## ---- 6. seed hygiene ---------------------------------------------------------
## Guards against (a) a seeded fit that is not reproducible and (b) the group
## draw perturbing the caller's RNG stream, which the reference script did
## (it called sample() directly).
fa <- mwperm_dyadic_het(dat$y, dat$d, x = dat$x, row = dat$i, col = dat$j,
                        n_flip = 5, n_reps = 3, seed = 11)
fb <- mwperm_dyadic_het(dat$y, dat$d, x = dat$x, row = dat$i, col = dat$j,
                        n_flip = 5, n_reps = 3, seed = 11)
stopifnot(isTRUE(same_fit(fa, fb)), length(fa$pvalues_rep) == 3L)
set.seed(99); before <- runif(3)
set.seed(99)
invisible(mwperm_dyadic_het(dat$y, dat$d, x = dat$x, row = dat$i, col = dat$j,
                            n_flip = 5, n_reps = 3, seed = 11))
invisible(build_flip_set(12, 12, 5, seed = 7))
stopifnot(identical(runif(3), before))
## different seeds draw different assignments (the test is a RANDOM test)
Fa <- build_flip_set(12, 12, 5, seed = 1)
Fb <- build_flip_set(12, 12, 5, seed = 2)
stopifnot(!identical(attr(Fa, "row_groups"), attr(Fb, "row_groups")))
## an unseeded call runs (and draws from the ambient stream)
set.seed(5)
stopifnot(inherits(build_flip_set(12, 12, 5), "list"))

## ---- 7. properties inherited from the engine --------------------------------
## Null-shift equivariance: the statistic is affine in b, so testing beta = b
## on y equals testing beta = 0 on y - D b, rep for rep (guards the claim that
## the confidence-set path carries over unchanged). And parallel == serial.
b0 <- 0.3
f_shift <- mwperm_dyadic_het(dat$y, dat$d, x = dat$x, row = dat$i, col = dat$j,
                             n_flip = 5, n_reps = 3, seed = 2, beta_null = b0,
                             conf_int = FALSE)
f_zero <- mwperm_dyadic_het(dat$y - b0 * dat$d, dat$d, x = dat$x, row = dat$i,
                            col = dat$j, n_flip = 5, n_reps = 3, seed = 2,
                            conf_int = FALSE)
stopifnot(identical(f_shift$pvalues_rep, f_zero$pvalues_rep))
if (.Platform$OS.type == "unix") {
  f_par <- mwperm_dyadic_het(dat$y, dat$d, x = dat$x, row = dat$i, col = dat$j,
                             n_flip = 5, n_reps = 3, seed = 11, n_cores = 2)
  stopifnot(isTRUE(same_fit(fa, f_par)))
}
## the methods all run on the object
out <- capture.output(print(f6))
stopifnot(any(grepl("Sign flips   : n_flip = 6 flip groups  (group order 2^5 = 32",
                    out, fixed = TRUE)),
          any(grepl("OLS estimate", out, fixed = TRUE)),
          any(grepl("IPT CI", out, fixed = TRUE)),
          !any(grepl("Permutations :", out, fixed = TRUE)))
ci6 <- confint(f6)
stopifnot(identical(dim(ci6), c(1L, 2L)),
          identical(attr(ci6, "method"), "IPT (inverted permutation test)"))
tf <- tempfile(fileext = ".pdf")
pdf(tf); plot(f6); dev.off()
stopifnot(file.exists(tf))
unlink(tf)

## ---- 8. the dispatcher: opt-in only, never auto-detected --------------------
## Guards against heteroskedasticity-robustness being "detected" (it cannot
## be: it leaves no trace in the clustering structure) and against the forced
## route returning anything other than the direct call.
chk <- mwperm_check(index = c("importer", "exporter"), data = trade_dyadic)
stopifnot(identical(chk$design, "dyadic"),
          any(grepl("dyadic_het", chk$alternatives, fixed = TRUE)),
          !any(grepl("dyadic_het", chk$notes, fixed = TRUE)))   # not a note
out_chk <- capture.output(print(chk))
stopifnot(any(grepl("design = \"dyadic_het\"", out_chk, fixed = TRUE)))
chk_h <- mwperm_check(index = c("importer", "exporter"), data = trade_dyadic,
                      design = "dyadic_het")
stopifnot(identical(chk_h$design, "dyadic_het"),
          identical(chk_h$n_flip_default, 6L),
          is.na(chk_h$K_default),
          identical(chk_h$p_floor, 1 / 32), isTRUE(chk_h$resolution_ok),
          identical(chk_h$levels_needed, 6L),
          grepl("mwperm_dyadic_het(", chk_h$call_str, fixed = TRUE))
out_h <- capture.output(print(chk_h))
stopifnot(any(grepl("default n_flip = 6", out_h, fixed = TRUE)),
          any(grepl("1/2^5 = 1/32", out_h, fixed = TRUE)))
## the diagnosis follows the same rule as the fit: alpha and aggregate move
## the default n_flip it reports, and the verdict stays attainable
chk_h2 <- mwperm_check(index = c("importer", "exporter"), data = trade_dyadic,
                       design = "dyadic_het", aggregate = "median2")
chk_h3 <- mwperm_check(index = c("importer", "exporter"), data = trade_dyadic,
                       design = "dyadic_het", alpha = 0.01)
stopifnot(identical(chk_h2$n_flip_default, 7L), isTRUE(chk_h2$resolution_ok),
          identical(chk_h2$p_floor, 2 / 64),
          identical(chk_h3$n_flip_default, 8L), isTRUE(chk_h3$resolution_ok),
          identical(chk_h3$p_floor, 1 / 128))
## a small design: the verdict names n_flip, not a cluster count
chk_s <- mwperm_check(index = list(i = g5$i, j = g5$j), design = "dyadic_het")
stopifnot(identical(chk_s$n_flip_default, 5L), isFALSE(chk_s$resolution_ok))
stopifnot(any(grepl("n_flip >= 6", capture.output(print(chk_s)), fixed = TRUE)))
## forced on an incomplete array: accepted (section 10), unlike "dyadic"
inc <- trade_dyadic[trade_dyadic$importer != trade_dyadic$exporter, ]
chk_i <- mwperm_check(index = c("importer", "exporter"), data = inc,
                      design = "dyadic_het")
stopifnot(identical(chk_i$design, "dyadic_het"),
          grepl("mwperm_dyadic_het(", chk_i$call_str, fixed = TRUE))
expect_err(mwperm_check(index = c("importer", "exporter"), data = inc,
                        design = "dyadic"), "complete")
## and the automatic route on an incomplete array (missing) still offers it
chk_a <- mwperm_check(index = c("importer", "exporter"), data = inc)
stopifnot(identical(chk_a$design, "missing"),
          any(grepl("dyadic_het", chk_a$alternatives, fixed = TRUE)))
## dispatch identity
f_dir <- with(trade_dyadic,
              mwperm_dyadic_het(y = log_trade, d = log_dist,
                                x = cbind(log_gdp_i, log_gdp_j),
                                row = importer, col = exporter,
                                n_flip = 6, n_reps = 2, seed = 1))
f_dis <- mwperm(y = "log_trade", d = "log_dist", x = c("log_gdp_i", "log_gdp_j"),
                index = c("importer", "exporter"), data = trade_dyadic,
                design = "dyadic_het", n_flip = 6, n_reps = 2, seed = 1,
                verbose = FALSE)
stopifnot(isTRUE(same_fit(f_dir, f_dis, skip = c("call", "auto"))),
          identical(f_dis$auto$design, "dyadic_het"))
f_frm <- mwperm_formula(log_trade ~ log_dist | log_gdp_i + log_gdp_j,
                        data = trade_dyadic, index = c("importer", "exporter"),
                        design = "dyadic_het", n_flip = 6, n_reps = 2, seed = 1,
                        verbose = FALSE)
stopifnot(isTRUE(same_fit(f_dir, f_frm, skip = c("call", "auto"))))
## `K` does not size this group, and `n_flip` does not apply elsewhere
expect_warn(mwperm(y = "log_trade", d = "log_dist", x = c("log_gdp_i", "log_gdp_j"),
                   index = c("importer", "exporter"), data = trade_dyadic,
                   design = "dyadic_het", K = 10, n_flip = 5, n_reps = 1,
                   seed = 1, conf_int = FALSE, verbose = FALSE),
            "`K` sizes a permutation group")
expect_warn(mwperm(y = "log_trade", d = "log_dist", x = c("log_gdp_i", "log_gdp_j"),
                   index = c("importer", "exporter"), data = trade_dyadic,
                   n_flip = 5, n_reps = 1, seed = 1, conf_int = FALSE,
                   verbose = FALSE),
            "`n_flip` applies to the dyadic_het design")

## ---- 9. the permutation designs are untouched --------------------------------
## Guards the numerics of every existing path through the shared engine: the
## seeded dyadic anchor (also pinned by test-dyadic.R and the golden
## baseline) and the permutation vocabulary of the resolution note.
fit_dy <- with(trade_dyadic,
               mwperm_dyadic(y = log_trade, d = log_dist,
                             x = cbind(log_gdp_i, log_gdp_j),
                             row = importer, col = exporter, seed = 1))
stopifnot(identical(fit_dy$pvalue, 0.025),
          abs(fit_dy$estimate - (-0.898503)) < 1e-6,
          abs(fit_dy$conf_int[1] - (-1.251351)) < 1e-5,
          abs(fit_dy$conf_int[2] - (-0.534468)) < 1e-5,
          is.null(fit_dy$n_flip),
          any(grepl("Permutations : K = 39", capture.output(print(fit_dy)),
                    fixed = TRUE)))
f_small <- mwperm_dyadic(y5, d5, row = g5$i, col = g5$j, seed = 1, n_reps = 1)
stopifnot(any(grepl(paste0("A 95% set needs K + 1 >= 20 -- that is, at least ",
                           "20 levels in the smallest permuted dimension."),
                    f_small$note, fixed = TRUE)))

## ---- 10. incomplete arrays ---------------------------------------------------
## A sign flip moves no observation, so the construction needs neither a
## complete array nor a biclique search: every observed cell is used and none
## is discarded. What it DOES need is that the kernel of the action stays
## {s, -s} on the observed cells. A flip group can be reachable only through
## cells that are all missing -- then its sign is free, the kernel doubles,
## and the 2^(n_flip - 1) representatives would contain duplicates while the
## fit reported a resolution it cannot attain. build_flip_set() therefore
## takes the observed cells and redraws until the row-group / column-group
## graph over them is connected.
##
## (a) on a complete array the observed-cell argument changes nothing: the
##     connectivity condition is implied by "every group used", so the draw
##     -- and every seeded result -- is bit-identical.
cells_all <- as.matrix(expand.grid(row = 1:6, col = 1:5))
for (sd in 1:5)
  stopifnot(identical(build_flip_set(6L, 5L, 4L, seed = sd),
                      build_flip_set(6L, 5L, 4L, seed = sd, cells = cells_all)))
## (b) a diagonal-only 4 x 4 array with n_flip = 4: without the guard most
##     assignments leave the group graph disconnected and the observation-
##     level sign vectors collapse onto fewer than 2^3 distinct ones.
cells_diag <- cbind(row = 1:4, col = 1:4)
for (sd in 1:40) {
  F <- build_flip_set(4L, 4L, 4L, seed = sd, cells = cells_diag)
  sv <- vapply(F, function(e) paste(e$row[cells_diag[, 1]] *
                                      e$col[cells_diag[, 2]], collapse = ""),
               "")
  stopifnot(length(F) == 8L, !anyDuplicated(sv))
}
expect_err(build_flip_set(4L, 4L, 4L, seed = 1, cells = cbind(1:2, 1:2)),
           "cells")                 # rows/cols outside 1..n_row / 1..n_col
## (c) the fit runs on an incomplete array, uses every observed cell, and is
##     Procedure 1 under the sign-flip group restricted to those cells: an
##     explicit orthonormal complement of [X | S_k X] built from scratch gives
##     the identical p-value.
set.seed(9)
inc2 <- inc[sample.int(nrow(inc), floor(0.8 * nrow(inc))), ]
f_inc <- with(inc2, mwperm_dyadic_het(y = log_trade, d = log_dist,
                                      x = cbind(log_gdp_i, log_gdp_j),
                                      row = importer, col = exporter,
                                      n_flip = 5, n_reps = 1, seed = 3,
                                      conf_int = FALSE))
stopifnot(f_inc$n_obs == nrow(inc2), f_inc$K == 15L,
          any(grepl("incomplete", f_inc$note, fixed = TRUE)),
          any(grepl("no cell is discarded", f_inc$note, fixed = TRUE)),
          any(grepl("independent of the errors given the covariates",
                    f_inc$note, fixed = TRUE)))
ri2 <- as.integer(factor(inc2$importer)); ci2 <- as.integer(factor(inc2$exporter))
F2 <- build_flip_set(max(ri2), max(ci2), 5L, seed = sub_seed(3L, 1L),
                     cells = cbind(ri2, ci2))
X2 <- cbind(1, inc2$log_gdp_i, inc2$log_gdp_j); D2 <- inc2$log_dist
y2 <- inc2$log_trade
null_basis <- function(M) {
  s <- svd(M, nu = nrow(M)); r <- sum(s$d > 1e-10 * s$d[1])
  s$u[, (r + 1):nrow(M), drop = FALSE]
}
a <- b <- numeric(length(F2) - 1L)
for (k in seq_along(a)) {
  sgn <- F2[[k + 1L]]$row[ri2] * F2[[k + 1L]]$col[ci2]
  V <- null_basis(cbind(X2, X2 * sgn)); P <- V %*% t(V)
  a[k] <- abs(sum(D2 * (P %*% y2))); b[k] <- abs(sum(D2 * (P %*% (y2 * sgn))))
}
stopifnot(identical(f_inc$pvalue, (1 + sum(b >= min(a))) / length(F2)))
## (d) dispatch identity on the incomplete array, and the complete-array
##     anchor is untouched (golden baseline; also section 9)
f_inc_dis <- mwperm(y = "log_trade", d = "log_dist",
                    x = c("log_gdp_i", "log_gdp_j"),
                    index = c("importer", "exporter"), data = inc2,
                    design = "dyadic_het", n_flip = 5, n_reps = 1, seed = 3,
                    conf_int = FALSE, verbose = FALSE)
stopifnot(isTRUE(same_fit(f_inc, f_inc_dis, skip = c("call", "auto"))))

passed("test-signflip.R")
