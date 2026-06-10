#' Downscale a coarse raster onto fine terrain by moving-window regression
#'
#' Fits a local linear regression of a coarse response on one or more coarse
#' predictors over a moving window, resamples the resulting intercept and slope
#' grids to the fine grid, and evaluates them on the fine-resolution predictors.
#' The output carries the fine-scale structure of the predictors with locally
#' varying coefficients.
#'
#' The relationship is given as a formula whose names refer to layers of `data`
#' (the coarse grid) and `onto` (the target grid). Because the response and the
#' coarse predictors are layers of a single `data` raster, they are guaranteed to
#' share a grid; the formula names match the predictors between `data` and `onto`,
#' so a missing or mis-ordered layer is an error rather than a silently wrong
#' result. Only bare layer names combined with `+` are supported, such as
#' `prec ~ elev + slope`. Transformations (`log(elev)`), interactions
#' (`elev:slope`), and `.` are rejected; create the derived layer first.
#'
#' In the common case there is one coarse response and a fine predictor such as a
#' digital elevation model, and no coarse predictor in hand. A predictor named in
#' the formula but absent from `data` is then derived by aggregating its `onto`
#' layer to the response grid with `aggregate`, so `topocast(prec ~ elev, data =
#' prec_coarse, onto = dem_fine, radius = 15)` works directly from a coarse
#' climate layer and a fine DEM.
#'
#' For a time series, supply `anomaly`: a stack of coarse periods. The baseline
#' relationship is fit once and each period's coarse anomaly, relative to
#' `baseline`, is carried onto the fine baseline. Use `type = "ratio"` for
#' non-negative variables such as precipitation and `type = "additive"` for
#' variables such as temperature.
#'
#' @param formula A two-sided formula of bare layer names, such as
#'   `prec ~ elev + slope`. The left-hand side names the coarse response layer in
#'   `data`; the right-hand side names the predictor layers.
#' @param data A `SpatRaster` on the coarse grid holding the response layer and,
#'   optionally, predictor layers named in `formula`. Any predictor not in `data`
#'   is derived from `onto`.
#' @param onto A `SpatRaster` on the target grid holding every predictor layer
#'   named in `formula`. Its grid defines the output.
#' @param radius Integer window radius in coarse cells; the window is a square of
#'   side `2 * radius + 1`.
#' @param aggregate Resampling method used to derive a coarse predictor from
#'   `onto` when it is not already a layer of `data`, passed to
#'   [terra::resample()]. Default `"average"`.
#' @param coefficients If `TRUE`, return the fitted layer together with the
#'   `(Intercept)` and per-predictor slope grids on the `onto` grid. Not supported
#'   with `anomaly`. Default `FALSE`.
#' @param anomaly Optional multi-layer `SpatRaster` on the coarse grid; each layer
#'   is one period to downscale relative to `baseline`. When supplied, the result
#'   has one layer per period.
#' @param baseline Optional single-layer `SpatRaster`, the coarse response baseline
#'   that `anomaly` is taken relative to. Defaults to the response layer of `data`.
#'   Ignored when `anomaly` is `NULL`.
#' @param type `"ratio"` (multiplicative) or `"additive"`; used only with
#'   `anomaly`.
#' @param method Resampling method for the coefficient grids, passed to
#'   [terra::resample()]. Default `"cubicspline"`.
#' @param min_cells,min_variance Passed to [window_regression()].
#'
#' @return A `SpatRaster` on the grid of `onto`. By default a single layer named
#'   for the response; one layer per period when `anomaly` is supplied; or the
#'   fitted layer plus `(Intercept)` and slope grids when `coefficients = TRUE`.
#'
#' @seealso [window_regression()] for the matrix engine.
#'
#' @examples
#' library(terra)
#' set.seed(1)
#' coarse <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20)
#' elevation <- setValues(coarse, runif(ncell(coarse), 0, 2000))
#' precipitation <- 800 - 0.1 * elevation
#' data <- c(precipitation, elevation)
#' names(data) <- c("prec", "elev")
#' terrain <- disagg(elevation, fact = 4, method = "bilinear")
#' names(terrain) <- "elev"
#'
#' # spatial downscale
#' fine <- topocast(prec ~ elev, data = data, onto = terrain, radius = 4)
#'
#' # one-DEM shortcut: the coarse predictor is derived from the fine DEM
#' fine2 <- topocast(prec ~ elev, data = data[["prec"]], onto = terrain, radius = 4)
#'
#' # return the local coefficient grids
#' coef_grids <- topocast(prec ~ elev, data = data, onto = terrain, radius = 4,
#'                        coefficients = TRUE)
#'
#' # time series: supply the periods
#' months <- precipitation * c(0.8, 1.2)
#' names(months) <- c("jan", "feb")
#' series <- topocast(prec ~ elev, data = data, onto = terrain, radius = 4,
#'                    anomaly = months, type = "ratio")
#'
#' @export
topocast <- function(formula, data, onto, radius,
                     aggregate = "average", coefficients = FALSE,
                     anomaly = NULL, baseline = NULL, type = c("ratio", "additive"),
                     method = "cubicspline", min_cells = 0L, min_variance = 1e-8) {
  type <- match.arg(type)
  if (coefficients && !is.null(anomaly))
    stop("`coefficients = TRUE` is not supported with `anomaly`; the coefficients ",
         "describe the baseline fit. Request them in a call without `anomaly`.")

  parsed <- parse_topo_formula(formula)
  if (!inherits(data, "SpatRaster")) stop("`data` must be a SpatRaster")
  if (!inherits(onto, "SpatRaster")) stop("`onto` must be a SpatRaster")
  onto <- harmonize_crs(data, onto)

  fit <- fit_windows(parsed, data, onto, radius = radius, aggregate = aggregate,
                     min_cells = min_cells, min_variance = min_variance)
  casted <- cast_onto(fit, onto, method = method)

  if (is.null(anomaly)) {
    if (coefficients) return(c(casted$fitted, casted$coefficients))
    return(casted$fitted)
  }

  if (!inherits(anomaly, "SpatRaster"))
    stop("`anomaly` must be a SpatRaster")
  if (is.null(baseline)) baseline <- data[[parsed$response]]
  carry_anomalies(casted$fitted, anomaly, baseline, type = type, method = method)
}

# --- internal: fit / cast / carry seam ------------------------------------

# Parse a `response ~ pred1 + pred2` formula into bare layer names, rejecting
# anything that is not a plain additive term.
parse_topo_formula <- function(formula) {
  if (!inherits(formula, "formula"))
    stop("`formula` must be a formula such as `prec ~ elev + slope`")

  # `.` means "all other columns" and would error in terms() without a data frame;
  # it has no meaning for named raster layers, so reject it up front.
  if ("." %in% all.names(formula[[length(formula)]]))
    stop(paste0("`formula` may only use bare layer names joined with `+`. ",
                "Transformations, interactions, and `.` are not supported; ",
                "create the derived layer before calling topocast()."))

  terms_obj <- stats::terms(formula)
  if (attr(terms_obj, "response") != 1L)
    stop("`formula` must have a response on the left-hand side, such as `prec ~ elev`")
  if (attr(terms_obj, "intercept") != 1L)
    stop("the intercept is always fit and cannot be removed from `formula`")

  response <- all.vars(formula[[2L]])
  if (length(response) != 1L)
    stop("the left-hand side of `formula` must be a single layer name")

  predictors <- attr(terms_obj, "term.labels")
  if (length(predictors) < 1L)
    stop("`formula` must name at least one predictor on the right-hand side")

  valid <- grepl("^[A-Za-z.][A-Za-z0-9._]*$", predictors) & predictors != "."
  if (!all(valid))
    stop(sprintf(
      paste0("`formula` may only use bare layer names joined with `+` (got: %s). ",
             "Transformations, interactions, and `.` are not supported; ",
             "create the derived layer before calling topocast()."),
      paste(predictors[!valid], collapse = ", ")))

  list(response = response, predictors = predictors)
}

# Treat two coordinate reference systems that share an EPSG code as equal even
# when their WKT strings differ, as happens with cross-source lon/lat data; fail
# with a message that names both systems otherwise.
harmonize_crs <- function(data, onto) {
  if (terra::same.crs(data, onto)) return(onto)

  desc_data <- terra::crs(data, describe = TRUE)
  desc_onto <- terra::crs(onto, describe = TRUE)
  same_code <- !is.na(desc_data$code) && !is.na(desc_onto$code) &&
    identical(desc_data$authority, desc_onto$authority) &&
    identical(desc_data$code, desc_onto$code)

  if (same_code) {
    onto <- terra::deepcopy(onto)
    terra::crs(onto) <- terra::crs(data)
    return(onto)
  }

  stop(sprintf(paste0(
    "`data` and `onto` do not share a coordinate reference system.\n",
    "  data: %s\n  onto: %s\n",
    "If these are the same projection with different WKT, align them with ",
    "`crs(onto) <- crs(data)` before calling topocast()."),
    crs_label(desc_data), crs_label(desc_onto)))
}

crs_label <- function(desc) {
  name <- if (!is.na(desc$name)) desc$name else "unknown CRS"
  if (!is.na(desc$code)) sprintf("%s (%s:%s)", name, desc$authority, desc$code) else name
}

# Pull named layers from a raster, with a clear error listing what is available.
select_layers <- function(raster, wanted, argument) {
  available <- names(raster)
  missing <- setdiff(wanted, available)
  if (length(missing))
    stop(sprintf(
      "`%s` is missing layer(s) named in `formula`: %s.\nAvailable layers: %s",
      argument, paste(missing, collapse = ", "), paste(available, collapse = ", ")))
  raster[[wanted]]
}

# Assemble the coarse predictors in formula order. A predictor that is a layer of
# `data` is used directly; one that is only in `onto` is aggregated to the
# response grid; one in neither is an error.
resolve_predictors <- function(predictors, data, onto, response_template, aggregate) {
  layers <- vector("list", length(predictors))
  for (i in seq_along(predictors)) {
    name <- predictors[i]
    if (name %in% names(data)) {
      layers[[i]] <- data[[name]]
    } else if (name %in% names(onto)) {
      layers[[i]] <- terra::resample(onto[[name]], response_template, method = aggregate)
    } else {
      stop(sprintf(paste0(
        "predictor `%s` named in `formula` is in neither `data` nor `onto`.\n",
        "  data layers: %s\n  onto layers: %s"),
        name, paste(names(data), collapse = ", "), paste(names(onto), collapse = ", ")))
    }
    names(layers[[i]]) <- name
  }
  terra::rast(layers)
}

# Fit the coarse coefficient grids. Returns the coefficients as a multi-layer
# SpatRaster (intercept plus one slope per predictor) on the coarse grid.
fit_windows <- function(parsed, data, onto, radius, aggregate, min_cells, min_variance) {
  response_raster <- select_layers(data, parsed$response, "data")
  predictor_raster <- resolve_predictors(parsed$predictors, data, onto,
                                         response_raster, aggregate)

  response_matrix <- raster_to_matrix(response_raster)
  predictor_matrices <- lapply(seq_len(terra::nlyr(predictor_raster)),
                               function(i) raster_to_matrix(predictor_raster[[i]]))
  regression <- window_regression(response_matrix, predictor_matrices,
                                  radius = radius, min_cells = min_cells,
                                  min_variance = min_variance)

  coefficients <- terra::rast(c(
    list(matrix_to_raster(regression$intercept, response_raster)),
    lapply(regression$slope, matrix_to_raster, template = response_raster)))
  names(coefficients) <- c("(Intercept)", parsed$predictors)
  list(coefficients = coefficients,
       response = parsed$response,
       predictors = parsed$predictors)
}

# Resample the coefficient grids to the target grid and evaluate them on the
# fine predictors: fitted = intercept + sum_j slope_j * predictor_j. Returns both
# the fitted layer and the resampled fine-grid coefficients.
cast_onto <- function(fit, onto, method) {
  fine_predictors <- select_layers(onto, fit$predictors, "onto")
  fine_coefficients <- terra::resample(fit$coefficients, fine_predictors[[1]], method = method)

  fitted <- fine_coefficients[["(Intercept)"]]
  for (predictor in fit$predictors)
    fitted <- fitted + fine_coefficients[[predictor]] * fine_predictors[[predictor]]
  names(fitted) <- fit$response
  list(fitted = fitted, coefficients = fine_coefficients)
}

# Carry each coarse period's anomaly, relative to a coarse baseline, onto the fine
# baseline. Ratio for non-negative variables, additive otherwise.
carry_anomalies <- function(fine_baseline, anomaly, baseline, type, method) {
  coarse_baseline <- terra::resample(baseline, fine_baseline, method = method)

  n_periods <- terra::nlyr(anomaly)
  periods <- vector("list", n_periods)
  for (period in seq_len(n_periods)) {
    coarse_value <- terra::resample(anomaly[[period]], fine_baseline, method = method)
    if (type == "ratio") {
      ratio <- coarse_value / coarse_baseline
      ratio <- terra::ifel(is.finite(ratio), ratio, NA)   # guard a zero or missing baseline
      periods[[period]] <- fine_baseline * ratio
    } else {
      periods[[period]] <- fine_baseline + (coarse_value - coarse_baseline)
    }
  }

  out <- terra::rast(periods)
  names(out) <- names(anomaly)
  out
}

# Raster to a wide matrix (row i, column j of the matrix is raster row i, col j).
raster_to_matrix <- function(raster) terra::as.matrix(raster, wide = TRUE)

# Wide matrix back to a raster on a template grid, in row-major cell order.
matrix_to_raster <- function(matrix, template) {
  terra::setValues(template, as.vector(t(matrix)))
}
