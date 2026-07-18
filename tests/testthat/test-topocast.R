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

test_that("the formula parses cbind() of bare names and rejects other left sides", {
  expect_equal(topocast:::parse_topo_formula(cbind(prec, temp) ~ elev + slope),
               list(response = c("prec", "temp"), predictors = c("elev", "slope")))
  expect_error(topocast:::parse_topo_formula(cbind(log(prec), temp) ~ elev),
               "bare layer names")
  expect_error(topocast:::parse_topo_formula(prec + temp ~ elev), "single bare layer name")
})

test_that("cbind() downscales several responses, each equal to its single call", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  temp <- 25 - 0.006 * grids$data[["elev"]]
  names(temp) <- "temp"
  data2 <- c(grids$data, temp)

  multi <- topocast(cbind(prec, temp) ~ elev + slope, data = data2,
                    onto = grids$terrain, radius = 3)
  expect_equal(names(multi), c("prec", "temp"))

  prec_one <- topocast(prec ~ elev + slope, data = data2, onto = grids$terrain, radius = 3)
  temp_one <- topocast(temp ~ elev + slope, data = data2, onto = grids$terrain, radius = 3)
  expect_equal(terra::values(multi[["prec"]]), terra::values(prec_one), tolerance = 1e-9)
  expect_equal(terra::values(multi[["temp"]]), terra::values(temp_one), tolerance = 1e-9)
})

test_that("several responses prefix their coefficient grids by response name", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  temp <- 25 - 0.006 * grids$data[["elev"]]; names(temp) <- "temp"
  out <- topocast(cbind(prec, temp) ~ elev, data = c(grids$data, temp),
                  onto = grids$terrain, radius = 3, coefficients = TRUE)
  expect_equal(names(out),
               c("prec", "prec.(Intercept)", "prec.elev",
                 "temp", "temp.(Intercept)", "temp.elev"))
})

test_that("diagnostics returns r.squared, residual.sd, and n.valid layers (issue #8)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  out <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                  radius = 3, diagnostics = TRUE)
  expect_equal(names(out), c("prec", "r.squared", "residual.sd", "n.valid"))

  r2 <- terra::values(out[["r.squared"]])
  r2 <- r2[is.finite(r2)]
  expect_true(all(r2 >= 0 & r2 <= 1))

  rsd <- terra::values(out[["residual.sd"]])
  rsd <- rsd[is.finite(rsd)]
  expect_true(all(rsd >= 0))

  nv <- terra::values(out[["n.valid"]])
  nv <- nv[is.finite(nv)]
  expect_true(all(nv > 0))
  expect_true(all(nv <= (2 * 3 + 1)^2))
})

test_that("with several responses, n.valid is shared (not prefixed) but r.squared/residual.sd are (issue #8)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  temp <- 25 - 0.006 * grids$data[["elev"]]; names(temp) <- "temp"
  out <- topocast(cbind(prec, temp) ~ elev, data = c(grids$data, temp),
                  onto = grids$terrain, radius = 3, diagnostics = TRUE)
  expect_equal(names(out),
               c("prec", "prec.r.squared", "prec.residual.sd",
                 "temp", "temp.r.squared", "temp.residual.sd", "n.valid"))
})

test_that("diagnostics adds residual.sd and n.valid alongside r.squared with anomaly (issue #8)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  out <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                  radius = 3, anomaly = grids$data[["prec"]], type = "ratio",
                  diagnostics = TRUE)
  expect_equal(names(out), c("prec", "r.squared", "residual.sd", "n.valid"))
})

test_that("clamp bounds the output to the coarse response range", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  rng <- range(terra::values(grids$data[["prec"]]), na.rm = TRUE)
  out <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                  radius = 3, clamp = TRUE)
  v <- terra::values(out); v <- v[is.finite(v)]
  expect_true(all(v >= rng[1] - 1e-9 & v <= rng[2] + 1e-9))
})

test_that("anomaly is rejected together with several responses", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  temp <- 25 - 0.006 * grids$data[["elev"]]; names(temp) <- "temp"
  expect_error(
    topocast(cbind(prec, temp) ~ elev, data = c(grids$data, temp), onto = grids$terrain,
             radius = 3, anomaly = grids$data[["prec"]]),
    "single-response")
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
