# topocast 0.0.5

* The moving-window fitting loop now runs on two threads while a package check is
  running, which R signals by setting `_R_CHECK_LIMIT_CORES_`, and on every
  available core otherwise. A check environment allows a package two cores, so the
  loop's previous fan-out over the whole machine spent more CPU time than the
  elapsed time a check budgets for it.

# topocast 0.0.4

* `topocast()` now brings every coarse grid a call fits -- each response's
  coefficients and diagnostics, plus the shared valid-cell count -- onto the target
  in one resample (grid target) or extract (point target), resolves the target
  predictors once, and assembles the coarse grids as a single multi-layer object
  rather than one terra object per coefficient. Each of those was previously done
  once per response, so a call repeated its whole coarse-to-target trip for every
  response it fit. The per-call fixed cost falls from about 0.58 s to about 0.11 s,
  a 25,000-point call with 19 responses from 1.28 s to 0.46 s, and the S4
  construction and dispatch that took over a fifth of a call no longer appears in
  its profile. With `anomaly`, the period stack likewise travels in one trip rather
  than one per period. Results are unchanged (#35).
* The moving-window engine's summed-area tables are now built and queried in
  double-double (extended) precision rather than plain `double`. A summed-area
  table entry near the far side of a large grid holds a sum over most of the
  grid's cells, so its magnitude grows with the grid size; on a plain-double
  table, a window's sum (extracted as the difference of four such entries) then
  loses precision proportional to that magnitude rather than to the window's own
  size, which by roughly a million cells could already turn a genuinely zero
  within-window variance into a small positive number that wrongly cleared
  `min_variance`, defeating the no-spread guard (#21).
* `arma::uword` is now explicitly 64-bit (`ARMA_64BIT_WORD`), and the moving
  window's row/column bounds are computed in a 64-bit type, closing two overflow
  classes left open after #13/#14 for a grid whose total cell count or single
  dimension is large enough to overflow a 32-bit computation (#31, #32).
* The C++ engine calls `Rcpp::checkUserInterrupt()` between chunks of rows (and
  during single-threaded setup), so a long-running fit on a large grid can now be
  interrupted rather than only checked for completion (#33).
* A response or predictor named `r.squared`, `residual.sd`, `n.valid`, or
  `(Intercept)` is now rejected at the formula-parsing stage; those names are
  reserved for `topocast()`'s own diagnostic and coefficient columns, and a
  response sharing one previously had its downscaled values silently overwritten
  by the diagnostic column of the same name when `diagnostics = TRUE` (#22).
* `n.valid` is now bounded to `[0, (2 * radius + 1)^2]` after resampling, the
  true range a window's valid-cell count can hold; like `r.squared`, it could
  previously be carried slightly outside its true range by some resampling
  methods (#24).
* `radius = 0` is now rejected: its one-cell window can never hold the `k + 1`
  valid cells a fit with `k` predictors needs, so it previously returned an
  all-`NA` result with no error or warning (#25).
* An invalid `aggregate` or grid-target `method` now errors with a message
  naming the actual `topocast()` argument, instead of terra's raw
  `'arg' should be one of ...` message, which named neither (#26).
* With `type = "ratio"`, a coarse cell where the baseline and the anomaly period
  are both exactly zero (e.g. a dry-season precipitation cell with no rain in
  either) is now treated as ratio `1`, "no change", rather than discarding the
  fine baseline the same way a true `x / 0` does (#27).
* A CRS mismatch error now falls back to a PROJ4 or WKT snippet when neither
  side carries an identifiable EPSG code (common for a custom or regional
  projection), instead of printing "unknown CRS" on both sides with no way to
  tell how they actually differ (#28).
* An entirely non-finite response/predictor stack now errors with a message
  naming the actual inputs (`data`/`onto` for `topocast()`, `y`/`x` for
  `window_regression()`) instead of the C++ engine's generic
  `"no finite cells to regress"` (#30, #34).
* `output = "raster"` no longer silently renames the `(Intercept)` coefficient
  column via the `raster` package's own `make.names()`-based `names<-`, which is
  now documented as an inherent limitation of that output class rather than left
  unmentioned (#23).
* `inst/CITATION` now reports the installed package version dynamically instead
  of a hardcoded string that had drifted out of date (#29).

* Several responses that share the predictors are downscaled in one call by naming
  them as `cbind(prec, tmin) ~ elev`. The moving-window design depends only on the
  predictors, so it is assembled and factored once and solved against every response;
  each extra response adds only a back-substitution. The result has one layer (or
  column) per response, and `coefficients`/`diagnostics` grids are prefixed by the
  response name. `window_regression()` likewise accepts a list of response matrices.
* `diagnostics = TRUE` returns an `r.squared` grid: the per-window coefficient of
  determination of the local fit, computed from the same summed-area sufficient
  statistics, mapping where the terrain relationship is strong. `window_regression()`
  now returns `r_squared` alongside `intercept` and `slope`.
* `diagnostics = TRUE` also returns `residual.sd` (the residual standard deviation of
  the local fit, in the response's own units) and `n.valid` (the count of valid
  coarse cells the window held), from the same sufficient statistics as `r.squared`.
  `n.valid` is shared across responses rather than prefixed, since the valid-cell mask
  is complete-case across them; `residual.sd` is prefixed like `r.squared` is.
  `window_regression()` likewise returns `residual_sd` and `n_valid` (#8).
* `clamp = TRUE` bounds the downscaled field to the observed range of the coarse
  response, a guard against the local linear fit extrapolating without limit where a
  fine predictor lies outside the range it was fit on.
* The per-cell window solve runs in parallel with OpenMP where the toolchain
  provides it; results are unchanged.
* `anomaly` and `baseline` accept a `Raster*` (raster) or `stars` object, like `data`
  and `onto` already did, instead of only a `SpatRaster` (#5).
* `radius` and `min_cells` are validated as non-negative whole numbers before
  reaching the C++ engine. A negative `radius` previously indexed the summed-area
  table out of bounds and crashed the R session rather than erroring (#4); the same
  validation now also rejects a value too large for a 32-bit integer, which
  previously overflowed to `NA` on coercion and crashed the session the same way (#9).
* `min_variance` is validated as a non-negative number. `NA` or a negative value
  previously compared as false against every window variance, silently disabling
  the documented guard against a predictor with no spread instead of being
  rejected (#10).
* `window_regression()` rejects a partially-named list of responses instead of
  silently making the unnamed ones' results unreachable by name (#6).
* `window_regression()` rejects a list of responses with a repeated name, which
  had the same unreachable-by-name problem as #6 but was not caught by it (#15).
* A `min_cells` close to `.Machine$integer.max` no longer overflows the C++
  engine's valid-cell threshold; previously the overflow silently disabled the
  guard instead of enforcing it, so windows with far fewer than the requested
  minimum could still be fit (#13).
* A `radius` close to `.Machine$integer.max` no longer crashes the R session;
  it is now clamped to the grid's span, which covers every cell in every
  window anyway and so changes no result (#14).
* `diagnostics = TRUE`'s `n.valid` (and `window_regression()`'s `n_valid`) now
  reports a window's valid-cell count whenever the window held enough cells to
  attempt a fit, rather than only when the fit also succeeded; a degenerate
  (no-spread predictor) or singular window previously reported `n.valid = NA`
  indistinguishably from a window that never had enough cells (#16).
* A window whose predictors are exactly collinear (rank-deficient design) now
  reliably returns `NA` instead of a spurious least-squares coefficient;
  Armadillo's default solver silently falls back to an approximate solution on
  a detected-singular system rather than failing, which defeated the
  documented singular-design guard (#17).

# topocast 0.0.3

* `data` and `onto` accept `Raster*` (raster) and `stars` objects in addition to
  `SpatRaster`, and the result is returned in the class of `onto`. The new
  `output` argument requests a specific class.
* `onto` may be an `sf` or `SpatVector` of points: the fitted relationship is
  evaluated at each point and returned as a prediction column, with the points
  carrying the fine predictor values as attributes. This makes downscaling to
  station or plot locations a single call. The `coefficients` and `anomaly`
  results are returned as columns in the same way.

# topocast 0.0.2

* `topocast()` derives a coarse predictor from `onto` when the predictor named in
  the formula is not a layer of `data`, aggregating it to the response grid with
  the new `aggregate` argument. The one-DEM case is now a single call from a
  coarse climate layer and a fine elevation model (#1).
* `coefficients = TRUE` returns the fitted layer together with the `(Intercept)`
  and per-predictor slope grids on the `onto` grid, so the local relationship,
  such as a precipitation lapse rate, can be mapped (#2).
* Coordinate reference systems that share an EPSG code are treated as equal even
  when their WKT strings differ, as happens with cross-source lon/lat data. A
  genuine mismatch now names both systems and suggests how to align them (#3).

# topocast 0.0.1

* First release.
* `topocast()` downscales a coarse raster onto fine terrain by moving-window
  regression. The relationship is a formula of layer names (`prec ~ elev +
  slope`); the response and coarse predictors live in one `data` raster and the
  fine predictors in `onto`. Coefficient grids are estimated with summed-area
  tables, so the cost is independent of the window radius.
* A time series is downscaled by passing a stack of coarse periods as `anomaly`:
  the baseline climatology is downscaled once and each period's anomaly is carried
  onto it, `type = "ratio"` or `type = "additive"`.
* `window_regression()` exposes the terra-free matrix engine.
