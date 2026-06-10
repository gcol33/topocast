# Build a small coarse data raster (response + predictors) and a fine target.
make_grids <- function(seed = 15) {
  set.seed(seed)
  coarse <- terra::rast(nrows = 12, ncols = 12, xmin = 0, xmax = 12, ymin = 0, ymax = 12)
  elevation <- terra::setValues(coarse, runif(terra::ncell(coarse), 0, 2000))
  slope <- terra::setValues(coarse, runif(terra::ncell(coarse), 0, 30))
  precipitation <- 800 - 0.1 * elevation + 2 * slope
  data <- c(precipitation, elevation, slope)
  names(data) <- c("prec", "elev", "slope")
  terrain <- terra::disagg(c(elevation, slope), fact = 4, method = "bilinear")
  names(terrain) <- c("elev", "slope")
  list(data = data, terrain = terrain)
}

test_that("topocast round-trips through terra onto the fine grid", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  out <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain, radius = 3)
  expect_s4_class(out, "SpatRaster")
  expect_equal(dim(out)[1:2], dim(grids$terrain)[1:2])
  expect_equal(names(out), "prec")
  expect_false(all(is.na(terra::values(out))))
})

test_that("the formula rejects non-additive and ill-formed terms", {
  expect_error(topocast:::parse_topo_formula(~ elev), "left-hand side")
  expect_error(topocast:::parse_topo_formula(prec ~ elev - 1), "intercept")
  expect_error(topocast:::parse_topo_formula(prec ~ log(elev)), "bare layer names")
  expect_error(topocast:::parse_topo_formula(prec ~ elev:slope), "bare layer names")
  expect_error(topocast:::parse_topo_formula(prec ~ .), "bare layer names")
  expect_equal(topocast:::parse_topo_formula(prec ~ elev + slope),
               list(response = "prec", predictors = c("elev", "slope")))
})

test_that("a layer named in the formula but absent is a clear error", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  expect_error(
    topocast(prec ~ elev + twi, data = grids$data, onto = grids$terrain, radius = 3),
    "missing layer")
  onto_without_slope <- grids$terrain[["elev"]]
  expect_error(
    topocast(prec ~ elev + slope, data = grids$data, onto = onto_without_slope, radius = 3),
    "missing layer")
})

test_that("an anomaly equal to the baseline returns the baseline (ratio and additive)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  baseline <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain, radius = 3)
  coarse_prec <- grids$data[["prec"]]

  # An anomaly equal to the baseline is the identity: ratio *1, additive +0.
  ratio <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                    radius = 3, anomaly = coarse_prec, type = "ratio")
  additive <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                       radius = 3, anomaly = coarse_prec, type = "additive")
  expect_equal(terra::values(ratio), terra::values(baseline), tolerance = 1e-6)
  expect_equal(terra::values(additive), terra::values(baseline), tolerance = 1e-6)
  expect_equal(terra::nlyr(ratio), terra::nlyr(coarse_prec))
})

test_that("a multi-period anomaly returns one layer per period, named", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  periods <- c(grids$data[["prec"]] * 0.5, grids$data[["prec"]] * 1.5)
  names(periods) <- c("dry", "wet")
  out <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                  radius = 3, anomaly = periods, type = "ratio")
  expect_equal(terra::nlyr(out), 2L)
  expect_equal(names(out), c("dry", "wet"))
})

test_that("the ratio path guards a zero baseline", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  zero_baseline <- terra::setValues(grids$data[["prec"]], 0)
  period <- grids$data[["prec"]]
  out <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                  radius = 3, anomaly = period, baseline = zero_baseline, type = "ratio")
  expect_true(all(is.na(terra::values(out))))
})
