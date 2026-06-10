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
