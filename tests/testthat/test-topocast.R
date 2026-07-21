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

test_that("cbind() with a duplicate response name errors at the formula parser (issue #19)", {
  expect_error(
    topocast:::parse_topo_formula(cbind(prec, prec) ~ elev),
    "names `prec` more than once")
})

test_that("a duplicate cbind() response name errors from topocast(), not window_regression() (issue #19)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  expect_error(
    topocast(cbind(prec, prec) ~ elev, data = grids$data, onto = grids$terrain, radius = 3),
    "names `prec` more than once")
})

test_that("anomaly and baseline are CRS-harmonized against data, like onto (issue #18)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  baseline_out <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                           radius = 3, anomaly = grids$data[["prec"]], type = "ratio")

  # Same EPSG code via an alternate WKT (cosmetic name change only): accepted,
  # and produces the same result.
  alt_anomaly <- grids$data[["prec"]]
  alt_wkt <- sub('"WGS 84 / UTM zone 32N"', '"WGS 84 / UTM zone 32N (alt)"',
                terra::crs(alt_anomaly), fixed = TRUE)
  terra::crs(alt_anomaly) <- alt_wkt
  out <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                  radius = 3, anomaly = alt_anomaly, type = "ratio")
  expect_equal(terra::values(out), terra::values(baseline_out), tolerance = 1e-9)

  # A genuine CRS mismatch on anomaly errors, the same way a mismatched onto does (issue #3).
  mismatched_anomaly <- grids$data[["prec"]]
  terra::crs(mismatched_anomaly) <- "EPSG:3857"
  expect_error(
    topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
             radius = 3, anomaly = mismatched_anomaly, type = "ratio"),
    "do not share a coordinate reference system")

  # A genuine CRS mismatch on a user-supplied baseline errors too.
  mismatched_baseline <- grids$data[["prec"]]
  terra::crs(mismatched_baseline) <- "EPSG:3857"
  expect_error(
    topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
             radius = 3, anomaly = grids$data[["prec"]], baseline = mismatched_baseline,
             type = "ratio"),
    "do not share a coordinate reference system")
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

test_that("a zero baseline with a zero period is treated as no change, not NA (issue #27)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  baseline_out <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain, radius = 3)
  zero_baseline <- terra::setValues(grids$data[["prec"]], 0)
  zero_period   <- terra::setValues(grids$data[["prec"]], 0)
  out <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                  radius = 3, anomaly = zero_period, baseline = zero_baseline, type = "ratio")
  expect_equal(terra::values(out), terra::values(baseline_out), tolerance = 1e-9)
})

test_that("a response or predictor name reserved for topocast()'s own output is rejected (issue #22)", {
  expect_error(topocast:::parse_topo_formula(`r.squared` ~ elev), "reserves")
  expect_error(topocast:::parse_topo_formula(`residual.sd` ~ elev), "reserves")
  expect_error(topocast:::parse_topo_formula(`n.valid` ~ elev), "reserves")
  expect_error(topocast:::parse_topo_formula(`(Intercept)` ~ elev), "reserves")
  expect_error(topocast:::parse_topo_formula(prec ~ r.squared), "reserves")
  expect_error(topocast:::parse_topo_formula(prec ~ n.valid), "reserves")
})

test_that("topocast() rejects a response named like a diagnostic column instead of silently overwriting it (issue #22)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  renamed <- grids$data
  names(renamed)[names(renamed) == "prec"] <- "r.squared"
  expect_error(
    topocast(`r.squared` ~ elev + slope, data = renamed, onto = grids$terrain,
             radius = 3, diagnostics = TRUE),
    "reserves")
})

test_that("every output class but raster preserves the (Intercept) coefficient column name (issue #23)", {
  skip_if_not_installed("terra")
  skip_if_not_installed("raster")
  skip_if_not_installed("stars")
  grids <- make_grids()
  args <- list(formula = prec ~ elev + slope, data = grids$data, onto = grids$terrain,
              radius = 3, coefficients = TRUE)
  expect_true("(Intercept)" %in% names(do.call(topocast, c(args, list(output = "terra")))))

  # stars stores a multi-layer raster as one 3D array with a "band" dimension,
  # not several named attributes; the per-layer names live in that dimension's
  # values rather than in names(), which instead reports the (single, and here
  # misleadingly first-layer-named) attribute.
  out_stars <- do.call(topocast, c(args, list(output = "stars")))
  expect_true("(Intercept)" %in% stars::st_get_dimension_values(out_stars, "band"))

  # raster::names() runs every name through make.names() unconditionally, in
  # both its setter and its getter; there is no way for a `raster` object to
  # report "(Intercept)" verbatim, so this output class alone gets the
  # make.names()-mangled name instead, as documented on `?topocast`.
  out_raster <- do.call(topocast, c(args, list(output = "raster")))
  expect_true("X.Intercept." %in% names(out_raster))
})

test_that("n.valid is bounded to [0, window cell count] after resampling (issue #24)", {
  expect_equal(topocast:::clamp_count(c(-5, 3, 100), 49), c(0, 3, 49))
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 2, ncols = 2, vals = c(-5, 3, 60, 49))
  out <- topocast:::clamp_count(r, 49)
  expect_equal(as.numeric(terra::values(out)), c(0, 3, 49, 49))
})

test_that("radius = 0 is rejected rather than silently returning an all-NA result (issue #25)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  expect_error(
    topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain, radius = 0),
    "at least 1")
})

test_that("an invalid aggregate or grid-target method names the topocast argument (issue #26)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  response_only <- grids$data[["prec"]]
  expect_error(
    topocast(prec ~ elev + slope, data = response_only, onto = grids$terrain, radius = 3,
             aggregate = "bogus_method"),
    "`aggregate")
  expect_error(
    topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain, radius = 3,
             method = "bogus_method"),
    "`method")
})

test_that("a CRS mismatch with no EPSG code on either side still shows something to compare (issue #28)", {
  skip_if_not_installed("terra")
  base  <- terra::rast(nrows = 4, ncols = 4,
    crs = "+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=500000 +datum=WGS84 +units=m +no_defs")
  other <- terra::rast(nrows = 4, ncols = 4,
    crs = "+proj=tmerc +lat_0=0 +lon_0=12 +k=0.9996 +x_0=500000 +datum=WGS84 +units=m +no_defs")
  expect_error(topocast:::harmonize_crs(base, other), "tmerc")
})

test_that("aggregate uses the requested resample method to derive a coarse predictor (issue #30)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  response_only <- grids$data[["prec"]]
  out_avg  <- topocast(prec ~ elev + slope, data = response_only, onto = grids$terrain,
                       radius = 3, aggregate = "average")
  out_near <- topocast(prec ~ elev + slope, data = response_only, onto = grids$terrain,
                       radius = 3, aggregate = "near")
  expect_false(isTRUE(all.equal(terra::values(out_avg), terra::values(out_near))))
})

test_that("cbind() downscales three or more responses correctly (issue #30)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  temp <- 25 - 0.006 * grids$data[["elev"]]; names(temp) <- "temp"
  wind <- 5 + 0.001 * grids$data[["elev"]]; names(wind) <- "wind"
  data3 <- c(grids$data, temp, wind)

  multi <- topocast(cbind(prec, temp, wind) ~ elev + slope, data = data3,
                    onto = grids$terrain, radius = 3)
  expect_equal(names(multi), c("prec", "temp", "wind"))

  wind_one <- topocast(wind ~ elev + slope, data = data3, onto = grids$terrain, radius = 3)
  expect_equal(terra::values(multi[["wind"]]), terra::values(wind_one), tolerance = 1e-9)
})

test_that("topocast() surfaces a clear error for an entirely non-finite response (issue #30/#34)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  all_na_prec <- terra::setValues(grids$data[["prec"]], NA_real_)
  names(all_na_prec) <- "prec"
  all_na_data <- c(all_na_prec, grids$data[["elev"]], grids$data[["slope"]])
  expect_error(
    topocast(prec ~ elev + slope, data = all_na_data, onto = grids$terrain, radius = 3),
    "nothing to fit")
})

test_that("a predictor present in both data and onto is fit from data's coarse layer (issue #30)", {
  skip_if_not_installed("terra")
  grids <- make_grids()
  perturbed_terrain <- grids$terrain
  perturbed_terrain[["elev"]] <- perturbed_terrain[["elev"]] + 5000

  coef_normal    <- topocast(prec ~ elev + slope, data = grids$data, onto = grids$terrain,
                             radius = 3, coefficients = TRUE)
  coef_perturbed <- topocast(prec ~ elev + slope, data = grids$data, onto = perturbed_terrain,
                             radius = 3, coefficients = TRUE)
  # the coarse fit (the coefficient grids) is unaffected by onto's elev values,
  # since a predictor already present in `data` is always fit from data's coarse
  # layer; only the final evaluation (not compared here) uses onto's elev.
  expect_equal(terra::values(coef_normal[["elev"]]), terra::values(coef_perturbed[["elev"]]),
               tolerance = 1e-9)
  expect_equal(terra::values(coef_normal[["(Intercept)"]]), terra::values(coef_perturbed[["(Intercept)"]]),
               tolerance = 1e-9)
})

test_that("topocast() passes threads down and the result does not depend on it (issue #36)", {
  set.seed(36)
  coarse <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40,
                        crs = "EPSG:32632")
  elevation <- terra::setValues(coarse, runif(terra::ncell(coarse), 0, 2000))
  precipitation <- 800 - 0.1 * elevation
  data <- c(precipitation, elevation)
  names(data) <- c("prec", "elev")
  terrain <- terra::disagg(elevation, fact = 2, method = "bilinear")
  names(terrain) <- "elev"

  serial <- topocast(prec ~ elev, data = data, onto = terrain, radius = 4, threads = 1)
  every  <- topocast(prec ~ elev, data = data, onto = terrain, radius = 4, threads = NULL)
  expect_equal(terra::values(serial), terra::values(every))

  expect_error(topocast(prec ~ elev, data = data, onto = terrain, radius = 4, threads = 0),
               "at least 1")
})
