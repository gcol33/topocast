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
#' For a time series, supply `anomaly`: a stack of coarse periods. The baseline
#' relationship is fit once and each period's coarse anomaly, relative to
#' `baseline`, is carried onto the fine baseline. Use `type = "ratio"` for
#' non-negative variables such as precipitation and `type = "additive"` for
#' variables such as temperature.
#'
#' @param formula A two-sided formula of bare layer names, such as
#'   `prec ~ elev + slope`. The left-hand side names the coarse response layer in
#'   `data`; the right-hand side names the predictor layers in both `data` and
#'   `onto`.
#' @param data A `SpatRaster` on the coarse grid holding the response layer and
#'   every predictor layer named in `formula`.
#' @param onto A `SpatRaster` on the target grid holding every predictor layer
#'   named in `formula`. Its grid defines the output.
#' @param radius Integer window radius in coarse cells; the window is a square of
#'   side `2 * radius + 1`.
#' @param anomaly Optional multi-layer `SpatRaster` on the coarse grid; each layer
#'   is one period to downscale relative to `baseline`. When supplied, the result
#'   has one layer per period.
#' @param baseline Optional single-layer `SpatRaster`, the coarse response baseline
#'   that `anomaly` is taken relative to. Defaults to the response layer of `data`.
#'   Ignored when `anomaly` is `NULL`.
#' @param type `"ratio"` (multiplicative) or `"additive"`; used only with
#'   `anomaly`.
#' @param method Resampling method for the coarse grids, passed to
#'   [terra::resample()]. Default `"cubicspline"`.
#' @param min_cells,min_variance Passed to [window_regression()].
#'
#' @return A `SpatRaster` on the grid of `onto`: a single layer named for the
#'   response, or one layer per period of `anomaly` when that is supplied.
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
#' # time series: same verb, supply the periods
#' months <- precipitation * c(0.8, 1.2)
#' names(months) <- c("jan", "feb")
#' series <- topocast(prec ~ elev, data = data, onto = terrain, radius = 4,
#'                    anomaly = months, type = "ratio")
#'
#' @export
topocast <- function(formula, data, onto, radius,
                     anomaly = NULL, baseline = NULL, type = c("ratio", "additive"),
                     method = "cubicspline", min_cells = 0L, min_variance = 1e-8) {
  type <- match.arg(type)

  fit <- fit_windows(formula, data = data, radius = radius,
                     min_cells = min_cells, min_variance = min_variance)
  fine_baseline <- cast_onto(fit, onto = onto, method = method)

  if (is.null(anomaly)) return(fine_baseline)

  if (!inherits(anomaly, "SpatRaster"))
    stop("`anomaly` must be a SpatRaster")
  if (is.null(baseline)) baseline <- data[[fit$response]]
  carry_anomalies(fine_baseline, anomaly, baseline, type = type, method = method)
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

# Fit the coarse coefficient grids. Returns the coefficients as a multi-layer
# SpatRaster (intercept plus one slope per predictor) on the coarse grid.
fit_windows <- function(formula, data, radius, min_cells, min_variance) {
  parsed <- parse_topo_formula(formula)
  if (!inherits(data, "SpatRaster"))
    stop("`data` must be a SpatRaster")

  response_raster <- select_layers(data, parsed$response, "data")
  predictor_raster <- select_layers(data, parsed$predictors, "data")

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
# fine predictors: fitted = intercept + sum_j slope_j * predictor_j.
cast_onto <- function(fit, onto, method) {
  if (!inherits(onto, "SpatRaster"))
    stop("`onto` must be a SpatRaster")
  if (!terra::same.crs(fit$coefficients, onto))
    stop("`data` and `onto` must share a coordinate reference system")

  fine_predictors <- select_layers(onto, fit$predictors, "onto")
  fine_coefficients <- terra::resample(fit$coefficients, fine_predictors[[1]], method = method)

  fitted <- fine_coefficients[["(Intercept)"]]
  for (predictor in fit$predictors)
    fitted <- fitted + fine_coefficients[[predictor]] * fine_predictors[[predictor]]
  names(fitted) <- fit$response
  fitted
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
