# mwperm 0.1.1

Publication-standard plotting (no change to any test, seed, or p-value).

* `plot.mwperm()` gains a `type` argument: `"coef"` (OLS estimate against the
  inverted-test confidence set -- an interval with end caps for one
  coefficient, a forest of the joint region's marginal extents for several),
  `"region"` (joint confidence region, two coefficients), `"stability"` (the
  Monte-Carlo p-value diagnostic, restyled), and `"all"`. **The default
  figure changed:** `type = "auto"` now draws the flagship `"coef"` figure
  whenever a confidence set is stored, falling back to `"stability"`
  otherwise (previously a p-value histogram/barplot was always drawn).
* All figures share one style layer: Okabe-Ito colourblind-safe palette with
  marks/line types differing as well as colour (grayscale-legible), no
  top/right spines, and a standard annotation block (design, cluster counts,
  N, resolution 1/(K+1), the null, the p-value, and the decision at alpha).
  Style elements can be overridden by name, e.g.
  `plot(fit, col_estimate = "black")`; unknown `...` arguments now warn and
  are ignored instead of crashing the underlying graphics calls.
* `"null"` and `"profile"` types are reserved (permutation-null and
  test-inversion p(b) figures); they need fit-time storage that this version
  does not retain and currently fall back to the default figure with a
  message.
* New `mwperm_save()` writes any figure at journal dimensions (single-column
  3.5 in / double-column 7 in, >= 300 dpi; pdf/png/tiff/jpeg).
* `grDevices` added to Imports (still base R only). Plotting never mutates
  global graphics state (`par()` is restored on exit).

# mwperm 0.1.0

Initial release.

* Finite-sample-valid invariant permutation tests (Guo, Toulis & Wang, 2026)
  for regression coefficients under multi-way clustering:
  `mwperm_dyadic()`, `mwperm_threeway()`, `mwperm_panel()` (arbitrary common
  time trend), `mwperm_layout()` (replicated two-way layouts, with `L0`
  balancing), and `mwperm_missing()` (incomplete arrays via fully observed
  bicliques, greedy or exact solver).
* Unified entry points: `mwperm()` auto-detects the design from the
  clustering structure and dispatches (with identical results to the direct
  call); `mwperm_check()` prints the diagnosis without running anything.
  Assumption-dependent forks (panel vs three-way; layout vs suppressed
  panel) default to the choice valid under the widest set of error processes
  and are announced with override instructions.
* Confidence sets by test inversion: an interval for one coefficient
  (with a guard that detects disconnected acceptance regions and widens to
  their hull with a note), a grid-based joint region for several.
* Opt-in parallelism via `n_cores` on all test functions (forked workers on
  Unix, PSOCK on Windows); results are identical to serial runs.
* S3 methods `print()`, `summary()`, `confint()`, `plot()`; synthetic
  example data `trade_dyadic` and `trade_panel` (generated reproducibly by
  `data-raw/make_data.R`).
