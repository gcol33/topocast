# Boundary adapters: accept raster, stars, sf, and SpatVector at the door, run the
# terra core unchanged, and return the class the user supplied (or asked for). The
# downscaling math lives in topocast.R; this file only translates classes.

# Require an optional package, with a message naming what it is needed for.
require_pkg <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("package '%s' is required for %s but is not installed.", pkg, what))
}

# The class token a gridded input's own class maps to, or NA for a non-grid input.
# The single place that knows which classes are grids, shared by as_grid() and
# as_target() so a class added to one is not forgotten in the other.
grid_class_of <- function(x) {
  if (inherits(x, "SpatRaster")) return("terra")
  if (inherits(x, c("RasterLayer", "RasterBrick", "RasterStack"))) return("raster")
  if (inherits(x, "stars")) return("stars")
  NA_character_
}

# Coerce a coarse gridded input to a SpatRaster. `data` is always a grid.
as_grid <- function(x, argument) {
  switch(grid_class_of(x),
    terra  = x,
    raster = { require_pkg("raster", "Raster* input"); terra::rast(x) },
    stars  = { require_pkg("stars", "stars input"); terra::rast(x) },
    stop(sprintf("`%s` must be a SpatRaster, Raster* (raster), or stars object.", argument)))
}

# Describe the target `onto` as either a grid (its geometry defines the output
# raster) or a set of point geometries (the fit is evaluated at each point). The
# original object and a class token are kept so the result can be returned in kind.
as_target <- function(onto) {
  class_token <- grid_class_of(onto)
  if (!is.na(class_token))
    return(list(kind = "grid", grid = as_grid(onto, "onto"), class = class_token))
  if (inherits(onto, "sf")) {
    require_pkg("sf", "sf input")
    vect <- terra::vect(onto)
    check_point_target(vect)
    return(list(kind = "vector", vect = vect, orig = onto, class = "sf"))
  }
  if (inherits(onto, "SpatVector")) {
    check_point_target(onto)
    return(list(kind = "vector", vect = onto, orig = onto, class = "spatvector"))
  }
  stop("`onto` must be a SpatRaster, Raster* (raster), stars, sf, or SpatVector object.")
}

# A vector target is evaluated point by point, so only point geometries are
# well-posed. Polygons need within-feature aggregation of the fine-grid result,
# which is a different operation; route them through a gridded `onto` instead.
check_point_target <- function(vect) {
  geom <- terra::geomtype(vect)
  if (geom != "points")
    stop(sprintf(paste0(
      "an sf/SpatVector `onto` must be points (got %s). topocast evaluates the local ",
      "regression at point locations; to summarise over polygons, downscale onto a ",
      "gridded `onto` and aggregate with terra::extract(result, polygons, fun = mean)."),
      geom))
}

# Harmonize the target's CRS against `data`, in place on whichever spatial object
# the target carries.
harmonize_target_crs <- function(data, target) {
  if (target$kind == "grid") target$grid <- harmonize_crs(data, target$grid)
  else                       target$vect <- harmonize_crs(data, target$vect)
  target
}

# Bring a coarse SpatRaster onto the target: resample to the fine grid, or
# interpolate at the point geometries. Returns a SpatRaster (grid) or a data frame
# whose columns are the coarse layers (points); both are indexed by layer name.
bring_onto_target <- function(coarse, target, method) {
  if (target$kind == "grid")
    return(terra::resample(coarse, target$grid[[1]], method = method))
  terra::extract(coarse, target$vect, method = "bilinear", ID = FALSE)
}

# Collapse a single-layer bring to a plain value: a one-layer SpatRaster (grid) or
# a numeric vector (points), so anomaly arithmetic is the same for both.
as_value <- function(x) if (inherits(x, "data.frame")) x[[1L]] else x

# The fine predictors named in the formula: layers of the grid target, or
# attribute columns of the point target. Both are indexed by name downstream.
target_predictors <- function(target, predictors) {
  if (target$kind == "grid")
    return(select_layers(target$grid, predictors, "onto"))
  values <- as.data.frame(target$vect)
  missing <- setdiff(predictors, names(values))
  if (length(missing))
    stop(sprintf(paste0(
      "`onto` (an sf/SpatVector target) is missing predictor attribute(s) named in ",
      "`formula`: %s.\nAvailable attributes: %s"),
      paste(missing, collapse = ", "), paste(names(values), collapse = ", ")))
  values
}

# fitted = intercept + sum_j slope_j * predictor_j. Works on SpatRaster layers and
# on numeric vectors alike, so the grid and point targets share one expression.
eval_fitted <- function(intercept, slopes, predictors) {
  fitted <- intercept
  for (name in names(slopes)) fitted <- fitted + slopes[[name]] * predictors[[name]]
  fitted
}

# A ratio with a zero or missing denominator set to NA, for SpatRaster or numeric.
safe_ratio <- function(numerator, denominator) {
  ratio <- numerator / denominator
  if (inherits(ratio, "SpatRaster")) return(terra::ifel(is.finite(ratio), ratio, NA))
  ifelse(is.finite(ratio), ratio, NA_real_)
}

# Clamp a value to a numeric range, for SpatRaster or numeric. A non-finite range
# (an all-NA response layer) is a no-op.
clamp_values <- function(x, range, target) {
  lo <- range[1L]; hi <- range[2L]
  if (!is.finite(lo) || !is.finite(hi)) return(x)
  if (target$kind == "grid") terra::clamp(x, lo, hi, values = TRUE)
  else pmin(pmax(x, lo), hi)
}

# Clamp to the unit interval, for SpatRaster or numeric. Resampling the coarse
# R-squared grid can carry it slightly outside [0, 1]; the fit value cannot.
clamp_unit <- function(x) {
  if (inherits(x, "SpatRaster")) return(terra::clamp(x, 0, 1, values = TRUE))
  pmin(pmax(x, 0), 1)
}

# Clamp to non-negative, for SpatRaster or numeric. Resampling the coarse residual-SD
# or valid-cell-count grid can carry it slightly below zero; neither value can be.
clamp_nonneg <- function(x) {
  if (inherits(x, "SpatRaster")) return(terra::clamp(x, lower = 0, values = TRUE))
  pmax(x, 0)
}

# The observed range of a coarse response layer, the bounds the downscaled field is
# clamped to when `clamp = TRUE`.
response_range <- function(data, response) {
  range(terra::values(data[[response]]), na.rm = TRUE)
}

# Prefix a named list of columns with `<response>.`, so several responses' coefficient
# or diagnostic grids do not collide in the output.
prefix_names <- function(cols, prefix) {
  stats::setNames(cols, paste0(prefix, ".", names(cols)))
}

# Assemble the named result columns into the target's representation and return it
# in the requested (or input-matching) class. `cols` are SpatRaster layers for a
# grid target and numeric vectors for a point target.
finalize <- function(target, cols, output) {
  class_out <- resolve_output_class(target, output)
  if (target$kind == "grid") {
    out <- terra::rast(cols)
    names(out) <- names(cols)
    return(as_grid_class(out, class_out))
  }
  frame <- as.data.frame(cols, check.names = FALSE)
  as_vector_class(target, frame, class_out)
}

# Pick the output class: follow the input by default, else the requested token,
# validated against what the target kind can produce.
resolve_output_class <- function(target, output) {
  if (is.null(output)) return(target$class)
  output <- match.arg(output, c("terra", "raster", "stars", "sf", "spatvector", "data.frame"))
  grid_ok   <- c("terra", "raster", "stars")
  vector_ok <- c("terra", "sf", "spatvector", "data.frame")
  allowed <- if (target$kind == "grid") grid_ok else vector_ok
  if (!output %in% allowed)
    stop(sprintf("`output = \"%s\"` is not available for a %s target; use one of: %s.",
                 output, if (target$kind == "grid") "raster/grid" else "sf/point",
                 paste(allowed, collapse = ", ")))
  output
}

# SpatRaster -> requested grid class.
as_grid_class <- function(x, class_out) {
  switch(class_out,
    terra  = x,
    raster = {
      require_pkg("raster", "Raster* output")
      if (terra::nlyr(x) == 1L) raster::raster(x) else raster::brick(x)
    },
    stars  = {
      require_pkg("stars", "stars output")
      stars::st_as_stars(x)
    },
    stop(sprintf("unsupported grid output class: %s", class_out)))
}

# Prediction columns -> requested vector class: appended to the input sf/SpatVector
# (geometry intact), or a plain data frame of coordinates plus predictions.
as_vector_class <- function(target, frame, class_out) {
  switch(class_out,
    sf = {
      require_pkg("sf", "sf output")
      out <- target$orig
      for (name in names(frame)) out[[name]] <- frame[[name]]
      out
    },
    terra = ,
    spatvector = {
      out <- terra::deepcopy(target$vect)
      for (name in names(frame)) out[[name]] <- frame[[name]]
      out
    },
    data.frame = {
      coords <- as.data.frame(terra::crds(target$vect, df = TRUE))
      cbind(coords, frame)
    },
    stop(sprintf("unsupported vector output class: %s", class_out)))
}
