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
