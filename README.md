# topocast

Moving-window regression downscaling of raster data.

`topocast` downscales a coarse raster (for example a climate variable) onto finer
terrain by fitting a local linear regression of the response on one or more
fine-resolution predictors (for example elevation and other terrain indices)
within a moving window. The regression coefficients are estimated for every cell
with summed-area tables, so the cost is independent of the window size, then
resampled to the target resolution and evaluated on the fine-resolution
predictors.

## Install

```r
# install.packages("pak")
pak::pak("gcol33/topocast")
```

## Use

The relationship is a formula whose names refer to layers of the coarse `data`
raster and the fine `onto` raster. The response and its predictors live in one
coarse raster, so they share a grid by construction; the formula names match the
predictors between the coarse and fine grids.

```r
library(topocast)
library(terra)

coarse  <- c(prec_1km, elev_1km)          # response + predictors, one grid
names(coarse) <- c("prec", "elev")
terrain <- elev_100m                       # target grid
names(terrain) <- "elev"

fine <- topocast(prec ~ elev, data = coarse, onto = terrain, radius = 15)
```

Add terrain predictors by naming more layers:

```r
coarse  <- c(prec_1km, elev_1km, twi_1km, slope_1km)
names(coarse)  <- c("prec", "elev", "twi", "slope")
terrain <- c(elev_100m, twi_100m, slope_100m)
names(terrain) <- c("elev", "twi", "slope")

fine <- topocast(prec ~ elev + twi + slope, data = coarse, onto = terrain, radius = 15)
```

For a time series, pass a stack of coarse periods as `anomaly`. The baseline
climatology is downscaled once and each period's coarse anomaly is carried onto it
(`type = "ratio"` for precipitation, `type = "additive"` for temperature):

```r
series <- topocast(
  prec ~ elev, data = coarse, onto = terrain, radius = 15,
  anomaly = prec_monthly_1km, type = "ratio"
)
```

The low-level matrix engine is exposed as `window_regression()`.
