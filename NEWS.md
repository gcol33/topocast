# topocast 0.0.4

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
  table out of bounds and crashed the R session rather than erroring (#4).
* `window_regression()` rejects a partially-named list of responses instead of
  silently making the unnamed ones' results unreachable by name (#6).

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
