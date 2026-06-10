# topocast: Moving-Window Regression Downscaling of Raster Data

Downscales coarse-resolution raster data to a finer grid by fitting
local linear regressions of a response, such as a climate variable, on
one or more fine-resolution predictors, such as elevation and other
terrain indices, within a moving window. Regression coefficients are
estimated for every cell using summed-area tables, so the cost is
independent of the window size, then resampled to the target resolution
and applied to the fine-resolution predictors. Multiplicative and
additive anomaly application downscale time series relative to a
baseline climatology, following the regression-on-elevation approach
used for high-resolution climate surfaces (Karger et al. 2017)
[doi:10.1038/sdata.2017.122](https://doi.org/10.1038/sdata.2017.122) .

## See also

Useful links:

- <https://gillescolling.com/topocast/>

- <https://github.com/gcol33/topocast>

- Report bugs at <https://github.com/gcol33/topocast/issues>

## Author

**Maintainer**: Gilles Colling <gilles.colling051@gmail.com>
([ORCID](https://orcid.org/0000-0003-3070-6066)) \[copyright holder\]

Authors:

- Gilles Colling <gilles.colling051@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-3070-6066)) \[copyright holder\]
