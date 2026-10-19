# Moving-window linear regression over raster matrices

Fits, for every cell, a linear regression of one or more responses on a
shared set of predictors over a square moving window, and returns
full-resolution intercept and slope grids together with a per-cell
coefficient of determination. The fit uses summed-area tables, so the
cost does not grow with the window size. Cells are excluded from a
window where any response or predictor is non-finite; a cell is returned
as `NA` when its window holds fewer valid cells than the model needs or
a predictor has no spread. If every cell of every response and predictor
is non-finite, there is nothing to regress and the call errors instead
of returning an all-`NA` result.

## Usage

``` r
window_regression(
  y,
  x,
  radius,
  min_cells = 0L,
  min_variance = 1e-08,
  threads = NULL
)
```

## Arguments

- y:

  Numeric matrix, the response on the coarse grid, or a list of numeric
  matrices for several responses sharing the predictors and grid.

- x:

  Numeric matrix, or a list of numeric matrices, the predictor(s) on the
  same grid as `y`.

- radius:

  Integer window radius in cells, at least `1`; the window is a square
  of side `2 * radius + 1`. A radius of `0` is rejected: its one-cell
  window can never hold the `k + 1` valid cells a fit needs.

- min_cells:

  Integer, additional valid cells required in a window beyond the
  `k + 1` model terms (`k` predictors plus the intercept). Default `0`.

- min_variance:

  Numeric, the minimum within-window variance a predictor must have for
  the cell to be fit. Default `1e-8`.

- threads:

  Integer, the number of threads the fitting loop may use, or `NULL`
  (the default) for every core available. A request above the number of
  cores the OpenMP runtime offers is clamped to that number, and every
  request is clamped to two while a package check is running. Set
  `threads = 1` when calling from inside a worker that is itself
  parallelized over targets, so the workers do not each spawn a full
  thread team; the fit is a small share of a call, so more workers at
  one thread each generally use the cores better than one worker at
  many. Without OpenMP the fit is single-threaded and this argument has
  no effect.

## Value

A list with `intercept` (a numeric matrix), `slope` (a list of numeric
matrices, one per predictor), `r_squared` (a numeric matrix),
`residual_sd` (a numeric matrix, the residual standard deviation of the
local fit; `NA` where the window has no residual degrees of freedom),
and `n_valid` (a numeric matrix, the count of valid cells the window
held), each the same size as `y`. `n_valid` is shared across responses,
since the valid-cell mask is complete-case across them. When `y` is a
list of responses, `intercept`, `r_squared`, and `residual_sd` are lists
of matrices and `slope` is a list of per-response slope lists, named for
the responses.

## Details

Several responses share the predictors, so the window design is
assembled and factored once and solved against every response, and each
extra response costs only a back-substitution. The valid-cell mask is
therefore complete-case across the responses: a cell enters a window
only where every response and predictor is finite.

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

# several responses sharing the predictor: the design is factored once
rain <- 800 - 0.1 * elevation
temp <- 15 - 0.006 * elevation
both <- window_regression(list(rain = rain, temp = temp), elevation, radius = 3)
names(both$slope)
```
