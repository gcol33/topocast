# Exact linear data so the window regression recovers a constant intercept (800)
# and slope (-0.1) everywhere; the point path then evaluates exactly and can be
# checked against the grid path.
make_io_grids <- function(seed = 21) {
  set.seed(seed)
  coarse <- terra::rast(nrows = 16, ncols = 16, xmin = 0, xmax = 16, ymin = 0, ymax = 16,
                        crs = "EPSG:32632")
  elevation <- terra::setValues(coarse, runif(terra::ncell(coarse), 0, 2000))
  precipitation <- 800 - 0.1 * elevation
  data <- c(precipitation, elevation)
  names(data) <- c("prec", "elev")
  terrain <- terra::disagg(elevation, fact = 4, method = "bilinear")
  names(terrain) <- "elev"
  list(data = data, terrain = terrain)
}

# Interior points carrying the bilinearly-interpolated elevation as an attribute.
make_io_points <- function(terrain) {
  coords <- cbind(x = c(4.3, 6.7, 8.1, 10.9, 12.5), y = c(5.2, 7.8, 9.4, 11.1, 6.3))
  elev <- terra::extract(terrain, coords, method = "bilinear")[["elev"]]
  sf::st_as_sf(data.frame(x = coords[, 1], y = coords[, 2], elev = elev),
               coords = c("x", "y"), crs = 32632)
}

test_that("a Raster* input round-trips to a Raster* output with matching values", {
  skip_if_not_installed("terra")
  skip_if_not_installed("raster")
  grids <- make_io_grids()
  terra_out <- topocast(prec ~ elev, data = grids$data, onto = grids$terrain, radius = 4)

  out <- topocast(prec ~ elev, data = raster::stack(grids$data),
                  onto = raster::stack(grids$terrain), radius = 4)
  expect_true(inherits(out, c("RasterLayer", "RasterStack", "RasterBrick")))
  expect_equal(as.numeric(terra::values(terra::rast(out))),
               as.numeric(terra::values(terra_out)), tolerance = 1e-6)
})

test_that("a stars input round-trips to a stars output with matching values", {
  skip_if_not_installed("terra")
  skip_if_not_installed("stars")
  grids <- make_io_grids()
  terra_out <- topocast(prec ~ elev, data = grids$data, onto = grids$terrain, radius = 4)

  out <- topocast(prec ~ elev, data = stars::st_as_stars(grids$data),
                  onto = stars::st_as_stars(grids$terrain), radius = 4)
  expect_s3_class(out, "stars")
  expect_equal(as.numeric(terra::values(terra::rast(out))),
               as.numeric(terra::values(terra_out)), tolerance = 1e-6)
})

test_that("output = requests an explicit class for a grid target", {
  skip_if_not_installed("terra")
  skip_if_not_installed("raster")
  skip_if_not_installed("stars")
  grids <- make_io_grids()
  as_raster <- topocast(prec ~ elev, data = grids$data, onto = grids$terrain,
                        radius = 4, output = "raster")
  as_stars <- topocast(prec ~ elev, data = grids$data, onto = grids$terrain,
                       radius = 4, output = "stars")
  expect_true(inherits(as_raster, c("RasterLayer", "RasterStack", "RasterBrick")))
  expect_s3_class(as_stars, "stars")
})

test_that("an sf point onto returns predictions evaluated at the points", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  grids <- make_io_grids()
  pts <- make_io_points(grids$terrain)

  out <- topocast(prec ~ elev, data = grids$data, onto = pts, radius = 4)
  expect_s3_class(out, "sf")
  expect_true("prec" %in% names(out))
  # exact linear data: prediction is intercept + slope * elev
  expect_equal(out$prec, 800 - 0.1 * pts$elev, tolerance = 1e-6)

  # the point path agrees with extracting the grid result at the same points
  terra_out <- topocast(prec ~ elev, data = grids$data, onto = grids$terrain, radius = 4)
  grid_at_pts <- terra::extract(terra_out, terra::vect(pts), method = "bilinear", ID = FALSE)[[1]]
  expect_equal(out$prec, grid_at_pts, tolerance = 1e-6)
})

test_that("coefficients = TRUE on a point target returns coefficient columns", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  grids <- make_io_grids()
  pts <- make_io_points(grids$terrain)

  out <- topocast(prec ~ elev, data = grids$data, onto = pts, radius = 4, coefficients = TRUE)
  expect_true(all(c("prec", "(Intercept)", "elev") %in% names(out)))
  expect_equal(unname(out[["(Intercept)"]]), rep(800, nrow(pts)), tolerance = 1e-4)
  expect_equal(unname(out[["elev"]]), rep(-0.1, nrow(pts)), tolerance = 1e-4)
})

test_that("an anomaly equal to the baseline is the identity on a point target", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  grids <- make_io_grids()
  pts <- make_io_points(grids$terrain)

  baseline <- topocast(prec ~ elev, data = grids$data, onto = pts, radius = 4)
  out <- topocast(prec ~ elev, data = grids$data, onto = pts, radius = 4,
                  anomaly = grids$data[["prec"]], type = "ratio")
  expect_equal(out$prec, baseline$prec, tolerance = 1e-6)
})

test_that("output = 'data.frame' returns coordinates plus predictions", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  grids <- make_io_grids()
  pts <- make_io_points(grids$terrain)

  out <- topocast(prec ~ elev, data = grids$data, onto = pts, radius = 4, output = "data.frame")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("x", "y", "prec") %in% names(out)))
})

test_that("a polygon onto is rejected with a clear message", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  grids <- make_io_grids()
  poly <- sf::st_sf(elev = 1000,
                    geometry = sf::st_sfc(sf::st_polygon(list(rbind(
                      c(2, 2), c(6, 2), c(6, 6), c(2, 6), c(2, 2)))), crs = 32632))
  expect_error(topocast(prec ~ elev, data = grids$data, onto = poly, radius = 4),
               "must be points")
})

test_that("an output class incompatible with the target kind errors", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  grids <- make_io_grids()
  pts <- make_io_points(grids$terrain)
  expect_error(
    topocast(prec ~ elev, data = grids$data, onto = grids$terrain, radius = 4, output = "sf"),
    "not available")
  expect_error(
    topocast(prec ~ elev, data = grids$data, onto = pts, radius = 4, output = "raster"),
    "not available")
})

test_that("a point onto missing a predictor attribute is a clear error", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  grids <- make_io_grids()
  pts <- make_io_points(grids$terrain)
  pts$elev <- NULL
  expect_error(topocast(prec ~ elev, data = grids$data, onto = pts, radius = 4),
               "predictor attribute")
})

test_that("a non-spatial onto is rejected", {
  skip_if_not_installed("terra")
  grids <- make_io_grids()
  expect_error(topocast(prec ~ elev, data = grids$data, onto = 1:10, radius = 4),
               "SpatRaster")
})
