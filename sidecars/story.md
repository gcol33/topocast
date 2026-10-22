# Story: topocast

## Claim
topocast downscales a coarse raster onto a fine grid of predictors by fitting a linear regression
in a moving window around every coarse cell, with summed-area tables so the cost does not grow
with the window.

## Audience
Ecologists and climate users who hold a coarse climate product (precipitation, temperature) and a
fine elevation model or terrain indices, and want the variable at the fine resolution. They know
terra rasters; they do not need to know summed-area tables.

## Relation
- to: none as a parent. The regression step is the regression-on-elevation approach described for
  CHELSA (Karger et al. 2017).
- phrase: "the regression step behind high-resolution climate surfaces such as CHELSA"
- same: a local linear regression of the coarse variable on terrain, applied to the fine terrain
- differs: fitted in every window rather than once, any number of named predictors, several
  responses per call
- status: not a successor of anything
- home: vignettes/how-it-works.Rmd
- context: none

## Arcs

### A1 Coarse grid to fine grid in one call
- question: "how do I get my 1 km climate layer onto my 100 m DEM?"
- answer: a formula, the coarse `data`, the fine `onto` and a `radius`; a predictor missing from
  `data` is aggregated from `onto`
- payoff: the reader downscales one layer with one call
- evidence: README usage chunk; getting-started.Rmd
- limits: `data` and `onto` share a coordinate reference system; the output does not reproduce the
  coarse input when averaged back, because the coarse residual is not carried
- home: vignettes/getting-started.Rmd
- also: README, DESCRIPTION

### A2 Several predictors, responses and periods
- question: "can I add slope and wetness, downscale temperature too, and carry a monthly series?"
- answer: named predictors on the right of the formula, `cbind()` responses on the left sharing one
  factored design, `anomaly` with `type = "ratio"` or `"additive"` for a series
- payoff: the reader downscales a set of variables and periods against one terrain
- evidence: getting-started.Rmd sections on predictors, responses and time series
- limits: `anomaly` takes one response; `coefficients = TRUE` is not available with `anomaly`
- home: vignettes/getting-started.Rmd
- also: README

### A3 Where the terrain explains the variable
- question: "is the downscaled detail real, and where?"
- answer: `coefficients = TRUE` maps the local intercept and slopes; `diagnostics = TRUE` maps
  `r.squared`, `residual.sd` and `n.valid`; `clamp = TRUE` bounds extrapolation
- payoff: the reader checks the fit before using the output
- evidence: getting-started.Rmd, how-it-works.Rmd "Reporting fit quality"
- limits: collinear predictors give unstable slopes in small windows
- home: vignettes/how-it-works.Rmd
- also: README, getting-started.Rmd

### A4 Cost independent of the window
- question: "will a large radius be slow?"
- answer: summed-area tables reduce each window fit to four lookups per sufficient statistic; the
  matrix engine `window_regression()` runs in C++ with OpenMP threads
- payoff: the reader picks the radius for the science, not the run time
- evidence: how-it-works.Rmd "Constant cost with summed-area tables"
- limits: none stated
- home: vignettes/how-it-works.Rmd
- also: README, DESCRIPTION

## Non-goals
- No climate data download.
- No residual correction that conserves the coarse values (the anomaly path carries each period's
  coarse anomaly instead).

## Documents
- DESCRIPTION: listing text; claim, A2 (anomaly), A4, the Karger et al. citation
- README / landing: entry point; claim, first example, one section per arc, routing
- Get started (getting-started.Rmd): tutorial; one landscape through every option
- how-it-works.Rmd: explanation; the method, A3 and A4 home

## Open
- No equivalence test against `lm()` on single windows is stated in the docs.
