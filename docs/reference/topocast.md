# Downscale a coarse raster onto fine terrain by moving-window regression

Fits a local linear regression of a coarse response on one or more
coarse predictors over a moving window, resamples the resulting
intercept and slope grids to the fine grid, and evaluates them on the
fine-resolution predictors. The output carries the fine-scale structure
of the predictors with locally varying coefficients.

## Usage

``` r
topocast(
  formula,
  data,
  onto,
  radius,
  aggregate = "average",
  coefficients = FALSE,
  anomaly = NULL,
  baseline = NULL,
  type = c("ratio", "additive"),
  method = "cubicspline",
  min_cells = 0L,
  min_variance = 1e-08
)
```

## Arguments

- formula:

  A two-sided formula of bare layer names, such as
  `prec ~ elev + slope`. The left-hand side names the coarse response
  layer in `data`; the right-hand side names the predictor layers.

- data:

  A `SpatRaster` on the coarse grid holding the response layer and,
  optionally, predictor layers named in `formula`. Any predictor not in
  `data` is derived from `onto`.

- onto:

  A `SpatRaster` on the target grid holding every predictor layer named
  in `formula`. Its grid defines the output.

- radius:

  Integer window radius in coarse cells; the window is a square of side
  `2 * radius + 1`.

- aggregate:

  Resampling method used to derive a coarse predictor from `onto` when
  it is not already a layer of `data`, passed to
  [`terra::resample()`](https://rspatial.github.io/terra/reference/resample.html).
  Default `"average"`.

- coefficients:

  If `TRUE`, return the fitted layer together with the `(Intercept)` and
  per-predictor slope grids on the `onto` grid. Not supported with
  `anomaly`. Default `FALSE`.

- anomaly:

  Optional multi-layer `SpatRaster` on the coarse grid; each layer is
  one period to downscale relative to `baseline`. When supplied, the
  result has one layer per period.

- baseline:

  Optional single-layer `SpatRaster`, the coarse response baseline that
  `anomaly` is taken relative to. Defaults to the response layer of
  `data`. Ignored when `anomaly` is `NULL`.

- type:

  `"ratio"` (multiplicative) or `"additive"`; used only with `anomaly`.

- method:

  Resampling method for the coefficient grids, passed to
  [`terra::resample()`](https://rspatial.github.io/terra/reference/resample.html).
  Default `"cubicspline"`.

- min_cells, min_variance:

  Passed to
  [`window_regression()`](https://gillescolling.com/topocast/reference/window_regression.md).

## Value

A `SpatRaster` on the grid of `onto`. By default a single layer named
for the response; one layer per period when `anomaly` is supplied; or
the fitted layer plus `(Intercept)` and slope grids when
`coefficients = TRUE`.

## Details

The relationship is given as a formula whose names refer to layers of
`data` (the coarse grid) and `onto` (the target grid). Because the
response and the coarse predictors are layers of a single `data` raster,
they are guaranteed to share a grid; the formula names match the
predictors between `data` and `onto`, so a missing or mis-ordered layer
is an error rather than a silently wrong result. Only bare layer names
combined with `+` are supported, such as `prec ~ elev + slope`.
Transformations (`log(elev)`), interactions (`elev:slope`), and `.` are
rejected; create the derived layer first.

In the common case there is one coarse response and a fine predictor
such as a digital elevation model, and no coarse predictor in hand. A
predictor named in the formula but absent from `data` is then derived by
aggregating its `onto` layer to the response grid with `aggregate`, so
`topocast(prec ~ elev, data = prec_coarse, onto = dem_fine, radius = 15)`
works directly from a coarse climate layer and a fine DEM.

For a time series, supply `anomaly`: a stack of coarse periods. The
baseline relationship is fit once and each period's coarse anomaly,
relative to `baseline`, is carried onto the fine baseline. Use
`type = "ratio"` for non-negative variables such as precipitation and
`type = "additive"` for variables such as temperature.

## See also

[`window_regression()`](https://gillescolling.com/topocast/reference/window_regression.md)
for the matrix engine.

## Examples

``` r
library(terra)
set.seed(1)
coarse <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
elevation <- setValues(coarse, runif(ncell(coarse), 0, 2000))
precipitation <- 800 - 0.1 * elevation
data <- c(precipitation, elevation)
names(data) <- c("prec", "elev")
terrain <- disagg(elevation, fact = 4, method = "bilinear")
names(terrain) <- "elev"

# spatial downscale
fine <- topocast(prec ~ elev, data = data, onto = terrain, radius = 4)

# one-DEM shortcut: the coarse predictor is derived from the fine DEM
fine2 <- topocast(prec ~ elev, data = data[["prec"]], onto = terrain, radius = 4)

# return the local coefficient grids
coef_grids <- topocast(prec ~ elev, data = data, onto = terrain, radius = 4,
                       coefficients = TRUE)

# time series: supply the periods
months <- precipitation * c(0.8, 1.2)
names(months) <- c("jan", "feb")
series <- topocast(prec ~ elev, data = data, onto = terrain, radius = 4,
                   anomaly = months, type = "ratio")
```
