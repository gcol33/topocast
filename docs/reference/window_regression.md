# Moving-window linear regression over raster matrices

Fits, for every cell, a linear regression of the response on one or more
predictors over a square moving window, and returns full-resolution
intercept and slope grids. The fit uses summed-area tables, so the cost
does not grow with the window size. Cells are excluded from a window
where the response or any predictor is non-finite; a cell is returned as
`NA` when its window holds fewer valid cells than the model needs or a
predictor has no spread.

## Usage

``` r
window_regression(y, x, radius, min_cells = 0L, min_variance = 1e-08)
```

## Arguments

- y:

  Numeric matrix, the response on the coarse grid.

- x:

  Numeric matrix, or a list of numeric matrices, the predictor(s) on the
  same grid as `y`.

- radius:

  Integer window radius in cells; the window is a square of side
  `2 * radius + 1`.

- min_cells:

  Integer, additional valid cells required in a window beyond the
  `k + 1` model terms (`k` predictors plus the intercept). Default `0`.

- min_variance:

  Numeric, the minimum within-window variance a predictor must have for
  the cell to be fit. Default `1e-8`.

## Value

A list with `intercept` (a numeric matrix) and `slope` (a list of
numeric matrices, one per predictor), each the same size as `y`.

## Details

This is the matrix engine behind
[`topocast()`](https://gillescolling.com/topocast/reference/topocast.md);
it works on plain numeric matrices and does not depend on terra.

## See also

[`topocast()`](https://gillescolling.com/topocast/reference/topocast.md)
for the terra workflow.

## Examples

``` r
set.seed(1)
elevation <- matrix(runif(100, 0, 1000), 10, 10)
climate <- 50 + 0.01 * elevation
fit <- window_regression(climate, elevation, radius = 3)
str(fit)
```
