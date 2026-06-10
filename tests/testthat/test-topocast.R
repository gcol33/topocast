# Build a small coarse data raster (response + predictors) and a fine target.
make_grids <- function(seed = 15) {
  set.seed(seed)
  coarse <- terra::rast(nrows = 12, ncols = 12, xmin = 0, xmax = 12, ymin = 0, ymax = 12,
                        crs = "EPSG:32632")
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

test_that("a predictor in neither data nor onto is a clear error", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  expect_error(
    topocast(prec ~ elev + twi, data = grids$data, onto = grids$terrain, radius = 3),
    "neither")
})

test_that("a predictor present in data but absent from onto errors on onto", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  onto_without_slope <- grids$terrain[["elev"]]
  expect_error(
    topocast(prec ~ elev + slope, data = grids$data, onto = onto_without_slope, radius = 3),
    "missing layer")
})

test_that("a coarse predictor is derived from onto when absent from data (issue #1)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  response_only <- grids$data[["prec"]]   # no predictor layers in data
  out <- topocast(prec ~ elev + slope, data = response_only, onto = grids$terrain, radius = 3)
  expect_s4_class(out, "SpatRaster")
  expect_equal(dim(out)[1:2], dim(grids$terrain)[1:2])
  expect_false(all(is.na(terra::values(out))))
})

test_that("coefficients = TRUE returns the fitted layer plus coefficient grids (issue #2)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  out <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                  radius = 3, coefficients = TRUE)
  expect_equal(names(out), c("prec", "(Intercept)", "elev", "slope"))
  fitted_only <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                          radius = 3)
  expect_equal(terra::values(out[["prec"]]), terra::values(fitted_only), tolerance = 1e-6)
})

test_that("coefficients = TRUE is rejected together with anomaly", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  expect_error(
    topocast(prec ~ elev, data = grids$data, onto = grids$terrain, radius = 3,
             coefficients = TRUE, anomaly = grids$data[["prec"]]),
    "not supported with")
})

test_that("harmonize_crs accepts equal codes and rejects different ones (issue #3)", {
  skip_if_not_installed("terra")
  base <- terra::rast(nrows = 4, ncols = 4, crs = "EPSG:4326")

  # Identical CRS: returned unchanged.
  expect_s4_class(topocast:::harmonize_crs(base, base), "SpatRaster")

  # Same EPSG code via an alternate WKT: accepted without error.
  alt <- base
  terra::crs(alt) <- paste0(
    'GEOGCRS["WGS 84",DATUM["World Geodetic System 1984",',
    'ELLIPSOID["WGS 84",6378137,298.257223563]],CS[ellipsoidal,2],',
    'AXIS["longitude",east],AXIS["latitude",north],',
    'ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",4326]]')
  expect_s4_class(topocast:::harmonize_crs(base, alt), "SpatRaster")

  # Different code: a clear error naming both systems.
  other <- base
  terra::crs(other) <- "EPSG:3857"
  expect_error(topocast:::harmonize_crs(base, other),
               "do not share a coordinate reference system")
})

test_that("a genuine CRS mismatch errors at the topocast() level (issue #3)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  onto_other <- grids$terrain
  terra::crs(onto_other) <- "EPSG:3857"
  expect_error(
    topocast(prec ~ elev + slope, data = grids$data, onto = onto_other, radius = 3),
    "do not share a coordinate reference system")
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
