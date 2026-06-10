# topocast

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Moving-window regression downscaling of raster data.

`topocast` rewrites a coarse raster at the resolution of a finer one by learning
how the coarse variable tracks the terrain, locally, and applying that
relationship to the fine-resolution predictors. A coarse precipitation grid plus
a digital elevation model becomes a precipitation field at the resolution of the
elevation model. The relationship is fit once per neighbourhood with a
moving-window regression whose cost does not grow with the window size.

## Quick start

```r
library(topocast)
library(terra)

coarse  <- c(prec_1km, elev_1km)      # response + predictors, one coarse grid
names(coarse) <- c("prec", "elev")
terrain <- elev_100m                   # the predictor at the target resolution
names(terrain) <- "elev"

fine <- topocast(prec ~ elev, data = coarse, onto = terrain, radius = 15)
```

## Statement of need

A climate variable often tracks elevation through a roughly linear relationship,
so a high-resolution elevation model can carry the variable to a finer grid
through a regression. This is the step behind high-resolution climate surfaces
such as CHELSA (Karger et al. 2017). Two things limit the usual implementations:
the relationship is fit globally or with a single predictor, and a moving-window
fit costs more as the window grows.

`topocast` fits the regression locally in a window around every cell, takes any
number of named terrain predictors, and uses summed-area tables so the per-cell
cost is four lookups regardless of the radius. The relationship is a formula of
layer names, which matches predictors between the coarse and fine grids and turns
a misnamed layer into an error rather than a silently wrong result.

## Features

- **Formula interface.** `prec ~ elev + slope + twi` names the response and the
  predictors; the names match layers between `data` and `onto`.
- **Any terrain predictors.** Elevation, slope, aspect, topographic wetness, or
  any aligned covariate, fit jointly.
- **Cost independent of window size.** Summed-area tables reduce each window fit
  to four lookups per sufficient statistic, so a radius of 30 costs the same as a
  radius of 3.
- **Time series in one call.** Pass a stack of coarse periods as `anomaly` to
  downscale the baseline once and carry each period onto it, ratio or additive.
- **A terra-free engine.** `window_regression()` exposes the matrix kernel for
  testing and for callers who hold their data as matrices.

## Installation

```r
# install.packages("pak")
pak::pak("gcol33/topocast")
```

## Usage

Several predictors are named in the formula and matched by name:

```r
coarse  <- c(prec_1km, elev_1km, twi_1km, slope_1km)
names(coarse)  <- c("prec", "elev", "twi", "slope")
terrain <- c(elev_100m, twi_100m, slope_100m)
names(terrain) <- c("elev", "twi", "slope")

fine <- topocast(prec ~ elev + twi + slope, data = coarse, onto = terrain,
                 radius = 15)
```

A time series shares one terrain relationship across periods. Pass the periods as
`anomaly` (`type = "ratio"` for precipitation, `type = "additive"` for
temperature):

```r
series <- topocast(prec ~ elev, data = coarse, onto = terrain, radius = 15,
                   anomaly = prec_monthly_1km, type = "ratio")
```

## Documentation

- [Getting started](https://gillescolling.com/topocast/articles/getting-started.html)
- [How moving-window downscaling works](https://gillescolling.com/topocast/articles/how-it-works.html)
- [Full reference](https://gillescolling.com/topocast/reference/)

## Support

> "Software is like sex: it's better when it's free." — Linus Torvalds

I'm a PhD student who builds R packages in my free time because I believe good
tools should be free and open. I started these projects for my own work and
figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to say
thanks. It helps with my coffee addiction.

[![Buy Me A Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License and citation

MIT. If you use `topocast` in published work, please cite it:

```bibtex
@software{topocast,
  author = {Colling, Gilles},
  title  = {topocast: Moving-Window Regression Downscaling of Raster Data},
  year   = {2026},
  url    = {https://github.com/gcol33/topocast}
}
```
