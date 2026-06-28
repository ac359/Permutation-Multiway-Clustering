## data-raw/make_data.R
## Reproducibly generate the two example datasets shipped with mwperm:
##   trade_dyadic  -- a synthetic bilateral-trade cross-section (complete array)
##   trade_panel   -- a synthetic bilateral-trade panel (complete (i, j, t) array)
##
## Both are built so that (a) a known covariate effect is present and
## recoverable, (b) the errors carry strong exporter/importer (and, for the
## panel, time) cluster components plus heavy tails, so that naive OLS standard
## errors are badly miscalibrated while the invariant permutation test stays
## valid, and (c) a deliberately irrelevant covariate ("placebo", true
## coefficient 0) is available for false-positive checks.
##
## Run from the package root:  Rscript data-raw/make_data.R

set.seed(20260601)

## a pool of real ISO3 country codes, just for flavour ----------------------
iso_pool <- c("USA","CHN","DEU","JPN","GBR","FRA","IND","ITA","BRA","CAN",
              "RUS","KOR","ESP","AUS","MEX","IDN","NLD","SAU","TUR","CHE",
              "POL","SWE","BEL","THA","ARG","AUT","NOR","ARE","ISR","ZAF",
              "DNK","SGP","EGY","PHL","CHL","FIN","COL","PRT","NZL","GRC")

## ===========================================================================
## 1. trade_dyadic : cross-section, n countries, complete n x n array
## ===========================================================================
n <- 40L
iso <- iso_pool[seq_len(n)]

## node (country) attributes
log_gdp   <- rnorm(n, mean = 12, sd = 1.4)         # economic mass of each country
lon       <- runif(n, -120, 120)                   # a pseudo-longitude for distance
lat       <- runif(n,  -40,  60)                    # a pseudo-latitude

grid <- expand.grid(i = seq_len(n), j = seq_len(n))  # includes internal trade i = j
i <- grid$i; j <- grid$j; N <- nrow(grid)

## dyad covariates
gcdist <- sqrt((lon[i] - lon[j])^2 + (lat[i] - lat[j])^2)
log_dist <- log(gcdist + 50)                        # +50 so that i = j is finite
border   <- as.integer(abs(lat[i] - lat[j]) < 12 & abs(lon[i] - lon[j]) < 18 & i != j)
placebo  <- rnorm(N)                                # irrelevant dyad covariate

## cluster-structured, heavy-tailed errors: exporter + importer + idiosyncratic
eta <- rnorm(n, sd = 0.8)                           # importer (row) effects
xi  <- rnorm(n, sd = 0.8)                           # exporter (col) effects
u   <- rt(N, df = 4) * 0.6                           # heavy-tailed idiosyncratic

## true structural coefficients
b0 <- 2.0; b_dist <- -1.0; b_gi <- 0.7; b_gj <- 0.7; b_border <- 0.4; b_placebo <- 0.0

log_trade <- b0 + b_dist * log_dist + b_gi * log_gdp[i] + b_gj * log_gdp[j] +
  b_border * border + b_placebo * placebo + eta[i] + xi[j] + u

trade_dyadic <- data.frame(
  importer  = factor(iso[i], levels = iso),
  exporter  = factor(iso[j], levels = iso),
  log_trade = round(log_trade, 4),
  log_dist  = round(log_dist, 4),
  log_gdp_i = round(log_gdp[i], 4),
  log_gdp_j = round(log_gdp[j], 4),
  border    = border,
  placebo   = round(placebo, 4),
  stringsAsFactors = FALSE
)
attr(trade_dyadic, "true_coef") <-
  c(`(Intercept)` = b0, log_dist = b_dist, log_gdp_i = b_gi,
    log_gdp_j = b_gj, border = b_border, placebo = b_placebo)

## ===========================================================================
## 2. trade_panel : ni countries x ni countries x nt years, complete array
## ===========================================================================
set.seed(20260602)
ni <- 22L; nt <- 6L
isop <- iso_pool[seq_len(ni)]
years <- 2015:(2015 + nt - 1L)

log_gdp0 <- rnorm(ni, mean = 12, sd = 1.3)
growth   <- runif(ni, 0.00, 0.06)                   # country-specific growth rate
lonp <- runif(ni, -120, 120); latp <- runif(ni, -40, 60)

pg <- expand.grid(i = seq_len(ni), j = seq_len(ni), t = seq_len(nt))
pi_ <- pg$i; pj <- pg$j; pt <- pg$t; Np <- nrow(pg)

## time-varying node GDP (a nuisance covariate)
log_gdp_i <- log_gdp0[pi_] + growth[pi_] * (pt - 1L)
log_gdp_j <- log_gdp0[pj] + growth[pj] * (pt - 1L)

## free-trade agreement: switches on for some pairs partway through the panel
set.seed(424242)
pair_id   <- pmin(pi_, pj) * 100L + pmax(pi_, pj)
fta_pairs <- unique(pair_id)
fta_start <- setNames(sample(c(rep(99L, length(fta_pairs) %/% 2),          # never
                               sample(seq_len(nt), length(fta_pairs) - length(fta_pairs) %/% 2,
                                      replace = TRUE))),                    # some year
                      as.character(fta_pairs))
fta <- as.integer(pt >= fta_start[as.character(pair_id)] & pi_ != pj)
placebo_p <- rnorm(Np)

## errors: importer + exporter cluster effects, an ARBITRARY (non-exchangeable)
## time trend zeta_t, and heavy-tailed idiosyncratic noise
etap  <- rnorm(ni, sd = 0.7)
xip   <- rnorm(ni, sd = 0.7)
zeta  <- c(0, 0.4, 0.25, 0.9, 0.7, 1.3)[seq_len(nt)]    # irregular over time
up    <- rt(Np, df = 4) * 0.5

a0 <- 1.5; a_fta <- 0.5; a_gi <- 0.6; a_gj <- 0.6; a_placebo <- 0.0
dist_p <- log(sqrt((lonp[pi_] - lonp[pj])^2 + (latp[pi_] - latp[pj])^2) + 50)

log_trade_p <- a0 + a_fta * fta + a_gi * log_gdp_i + a_gj * log_gdp_j -
  0.8 * dist_p + a_placebo * placebo_p + etap[pi_] + xip[pj] + zeta[pt] + up

trade_panel <- data.frame(
  importer  = factor(isop[pi_], levels = isop),
  exporter  = factor(isop[pj], levels = isop),
  year      = years[pt],
  log_trade = round(log_trade_p, 4),
  fta       = fta,
  log_gdp_i = round(log_gdp_i, 4),
  log_gdp_j = round(log_gdp_j, 4),
  placebo   = round(placebo_p, 4),
  stringsAsFactors = FALSE
)
attr(trade_panel, "true_coef") <-
  c(`(Intercept)` = a0, fta = a_fta, log_gdp_i = a_gi, log_gdp_j = a_gj,
    placebo = a_placebo)

## ===========================================================================
## save to data/
## ===========================================================================
dir.create("data", showWarnings = FALSE)
save(trade_dyadic, file = "data/trade_dyadic.rda", version = 2, compress = "xz")
save(trade_panel,  file = "data/trade_panel.rda",  version = 2, compress = "xz")

cat(sprintf("trade_dyadic: %d rows (%d countries, complete %dx%d)\n",
            nrow(trade_dyadic), n, n, n))
cat(sprintf("trade_panel : %d rows (%d countries x %d years)\n",
            nrow(trade_panel), ni, nt))
