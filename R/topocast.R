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
#' Several responses that share the predictors are downscaled together by naming
#' them as `cbind(prec, tmin) ~ elev`. The moving-window design depends only on the
#' predictors, so it is fit once and solved against every response; the result has
#' one layer (or column) per response. With `coefficients` or `diagnostics`, the
#' extra grids are prefixed by the response name.
#'
#' For a time series, supply `anomaly`: a stack of coarse periods. The baseline
#' relationship is fit once and each period's coarse anomaly, relative to
#' `baseline`, is carried onto the fine baseline. Use `type = "ratio"` for
#' non-negative variables such as precipitation and `type = "additive"` for
#' variables such as temperature.
#'
#' `data` and `onto` may be any of the common spatial classes. A gridded input is
#' accepted as a `SpatRaster`, a `Raster*` object (raster), or a `stars` object.
#' The target `onto` may instead be a set of points as an `sf` or `SpatVector`
#' object, in which case the fitted relationship is evaluated at each point and a
#' prediction column is returned; the points must carry the fine predictor values
#' as attributes. By default the result is returned in the class of `onto`; set
#' `output` to request another.
#'
#' @param formula A two-sided formula of bare layer names, such as
#'   `prec ~ elev + slope`. The left-hand side names the coarse response layer in
#'   `data`, or `cbind(prec, tmin)` for several responses sharing the predictors;
#'   the right-hand side names the predictor layers.
#' @param data A gridded coarse input holding the response layer and, optionally,
#'   predictor layers named in `formula`: a `SpatRaster`, a `Raster*` (raster), or
#'   a `stars` object. Any predictor not in `data` is derived from `onto`.
#' @param onto The target. A gridded `SpatRaster`, `Raster*`, or `stars` object
#'   whose grid defines the output, holding every predictor layer named in
#'   `formula`; or an `sf`/`SpatVector` of points carrying those predictors as
#'   attributes, in which case the fit is evaluated at the points. With a point
#'   `onto` every predictor must be a layer of `data` (the derive-from-`onto`
#'   shortcut needs a grid).
#' @param radius Integer window radius in coarse cells; the window is a square of
#'   side `2 * radius + 1`.
#' @param aggregate Resampling method used to derive a coarse predictor from
#'   `onto` when it is not already a layer of `data`, passed to
#'   [terra::resample()]. Default `"average"`.
#' @param coefficients If `TRUE`, return the fitted layer together with the
#'   `(Intercept)` and per-predictor slope grids on the `onto` grid. Not supported
#'   with `anomaly`. Default `FALSE`.
#' @param diagnostics If `TRUE`, also return three grids (or columns) describing the
#'   local fit, each brought onto `onto`: `r.squared`, the per-window coefficient of
#'   determination, mapping where the terrain relationship is strong and where the
#'   downscaled field rests mostly on the coarse input; `residual.sd`, the residual
#'   standard deviation of the local fit in the response's own units; and `n.valid`,
#'   the count of valid coarse cells the window held, which tells a fit that is weak
#'   because the window barely cleared the minimum valid-cell requirement apart from
#'   one that is weak because the terrain relationship is genuinely noisy there.
#'   `n.valid` is reported once per call, not once per response, since the
#'   valid-cell mask is complete-case across every response fit together. Default
#'   `FALSE`.
#' @param anomaly Optional multi-layer coarse grid (a `SpatRaster`, `Raster*`, or
#'   `stars` object, like `data`); each layer is one period to downscale relative
#'   to `baseline`. When supplied, the result has one layer (or column) per period.
#' @param baseline Optional single-layer coarse grid (a `SpatRaster`, `Raster*`, or
#'   `stars` object, like `data`), the coarse response baseline that `anomaly` is
#'   taken relative to. Defaults to the response layer of `data`. Ignored when
#'   `anomaly` is `NULL`.
#' @param type `"ratio"` (multiplicative) or `"additive"`; used only with
#'   `anomaly`.
#' @param method Interpolation method for bringing the coefficient grids onto
#'   `onto`. For a grid `onto`, passed to [terra::resample()]; default
#'   `"cubicspline"`. For a point (sf/SpatVector) `onto`, passed to
#'   [terra::extract()], which only supports `"simple"` and `"bilinear"`;
#'   default `"bilinear"`. Default `NULL` uses the kind-appropriate default;
#'   an explicit value that `terra::extract()` does not support is an error
#'   for a point target.
#' @param output Optional output class, one of `"terra"`, `"raster"`, `"stars"`
#'   (grid targets) or `"terra"`, `"sf"`, `"spatvector"`, `"data.frame"` (point
#'   targets). Default `NULL` returns the result in the class of `onto`.
#' @param clamp If `TRUE`, bound the downscaled field to the observed range of the
#'   coarse response, a guard against the local linear fit extrapolating without
#'   limit where a fine predictor lies outside the range it was fit on. Default
#'   `FALSE`.
#' @param min_cells,min_variance Passed to [window_regression()].
#'
#' @return The downscaled result on the geometry of `onto`, in the class of `onto`
#'   or the class named by `output`. For a grid target: a single layer named for
#'   the response; one layer per response with a `cbind()` left-hand side; one layer
#'   per period when `anomaly` is supplied; the fitted layer plus `(Intercept)` and
#'   slope grids when `coefficients = TRUE`; and `r.squared` and `residual.sd` grids
#'   plus a shared `n.valid` grid when `diagnostics = TRUE`. With several responses
#'   the coefficient, `r.squared`, and `residual.sd` grids are prefixed by the
#'   response name; `n.valid` is not, since it is the same valid-cell mask for every
#'   response fit in the call. For a point target the same quantities are returned
#'   as prediction columns.
#'
#' @seealso [window_regression()] for the matrix engine.
#'
#' @examples
#' library(terra)
#' set.seed(1)
#' coarse <- rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 0, ymax = 20,
#'                crs = "EPSG:32632")
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
#' # several responses sharing the predictor, downscaled in one call
#' temperature <- 25 - 0.006 * elevation
#' data2 <- c(precipitation, temperature, elevation)
#' names(data2) <- c("prec", "temp", "elev")
#' both <- topocast(cbind(prec, temp) ~ elev, data = data2, onto = terrain, radius = 4)
#'
#' # diagnostics: the local fit quality as an r.squared grid
#' with_r2 <- topocast(prec ~ elev, data = data, onto = terrain, radius = 4,
#'                     diagnostics = TRUE)
#'
#' # time series: supply the periods
#' months <- precipitation * c(0.8, 1.2)
#' names(months) <- c("jan", "feb")
#' series <- topocast(prec ~ elev, data = data, onto = terrain, radius = 4,
#'                    anomaly = months, type = "ratio")
#'
#' # predict at point locations: onto is sf points carrying the predictor
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   plots <- sf::st_as_sf(
#'     data.frame(x = c(5, 10, 15), y = c(5, 10, 15), elev = c(500, 1000, 1500)),
#'     coords = c("x", "y"), crs = "EPSG:32632")
#'   at_plots <- topocast(prec ~ elev, data = data, onto = plots, radius = 4)
#' }
#'
#' @export
topocast <- function(formula, data, onto, radius,
                     aggregate = "average", coefficients = FALSE, diagnostics = FALSE,
                     anomaly = NULL, baseline = NULL, type = c("ratio", "additive"),
                     method = NULL, output = NULL, clamp = FALSE,
                     min_cells = 0L, min_variance = 1e-8) {
  type <- match.arg(type)
  parsed <- parse_topo_formula(formula)
  multi  <- length(parsed$response) > 1L

  if (coefficients && !is.null(anomaly))
    stop("`coefficients = TRUE` is not supported with `anomaly`; the coefficients ",
         "describe the baseline fit. Request them in a call without `anomaly`.")
  if (multi && !is.null(anomaly))
    stop("`anomaly` downscaling is single-response: the period stack is taken ",
         "relative to one baseline. Call topocast() once per response.")

  data   <- as_grid(data, "data")
  target <- as_target(onto)
  target <- harmonize_target_crs(data, target)
  method <- resolve_cast_method(method, target)

  onto_grid <- if (target$kind == "grid") target$grid else NULL
  fit <- fit_windows(parsed, data, onto_grid, radius = radius, aggregate = aggregate,
                     min_cells = min_cells, min_variance = min_variance)

  if (is.null(anomaly)) {
    cols <- list()
    for (rfit in fit$responses) {
      resp   <- rfit$response
      casted <- cast_onto(rfit, fit$predictors, target, method = method)
      fitted <- casted$fitted
      if (clamp) fitted <- clamp_values(fitted, response_range(data, resp), target)
      cols[[resp]] <- fitted
      if (coefficients)
        cols <- c(cols, if (multi) prefix_names(casted$coef_cols, resp) else casted$coef_cols)
      if (diagnostics) {
        cols[[if (multi) paste0(resp, ".r.squared")   else "r.squared"]]   <- casted$r_squared
        cols[[if (multi) paste0(resp, ".residual.sd") else "residual.sd"]] <- casted$residual_sd
      }
    }
    if (diagnostics) cols[["n.valid"]] <- cast_n_valid(fit$n_valid, target, method)
    return(finalize(target, cols, output))
  }

  anomaly <- harmonize_crs(data, as_grid(anomaly, "anomaly"))
  rfit   <- fit$responses[[1L]]
  resp   <- rfit$response
  casted <- cast_onto(rfit, fit$predictors, target, method = method)
  fine_baseline <- casted$fitted
  if (clamp) fine_baseline <- clamp_values(fine_baseline, response_range(data, resp), target)
  baseline <- if (is.null(baseline)) data[[resp]] else harmonize_crs(data, as_grid(baseline, "baseline"))
  cols <- carry_anomalies(fine_baseline, anomaly, baseline, target,
                          type = type, method = method)
  if (diagnostics) {
    cols[["r.squared"]]   <- casted$r_squared
    cols[["residual.sd"]] <- casted$residual_sd
    cols[["n.valid"]]     <- cast_n_valid(fit$n_valid, target, method)
  }
  finalize(target, cols, output)
}

# --- internal: fit / cast / carry seam ------------------------------------

# Parse a `response ~ pred1 + pred2` formula into bare layer names, rejecting
# anything that is not a plain additive term. The left-hand side is a single bare
# name, or `cbind(r1, r2, ...)` of bare names for several responses sharing the
# predictors.
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

  response <- parse_response(formula[[2L]])

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

# Parse the left-hand side: a single bare name, or `cbind(r1, r2, ...)` of bare
# names for several responses. Anything else (a transformation, a literal) is an
# error, mirroring the right-hand-side rule.
parse_response <- function(lhs) {
  if (is.name(lhs)) return(as.character(lhs))
  if (is.call(lhs) && identical(lhs[[1L]], as.name("cbind"))) {
    args <- as.list(lhs)[-1L]
    if (length(args) < 1L)
      stop("`cbind()` on the left-hand side must name at least one response layer")
    if (!all(vapply(args, is.name, logical(1))))
      stop(paste0("the left-hand side may only use bare layer names; use ",
                  "`cbind(r1, r2)` of bare names for several responses, such as ",
                  "`cbind(prec, tmin) ~ elev`."))
    names <- vapply(args, as.character, character(1))
    dup <- names[duplicated(names)]
    if (length(dup))
      stop(sprintf(paste0(
        "the left-hand side of `cbind()` names `%s` more than once; ",
        "each response must have a unique name."), dup[1L]))
    return(names)
  }
  stop(paste0("the left-hand side of `formula` must be a single bare layer name, ",
              "or `cbind(...)` of bare layer names for several responses."))
}

# Treat two coordinate reference systems that share an EPSG code as equal even
# when their WKT strings differ, as happens with cross-source lon/lat data; fail
# with a message that names both systems otherwise. Works on SpatRaster and
# SpatVector targets alike.
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
# `data` is used directly; one that is only in the `onto` grid is aggregated to the
# response grid. With a point target there is no `onto` grid, so a predictor must
# be a layer of `data`; anything missing is an error.
resolve_predictors <- function(predictors, data, onto_grid, response_template, aggregate) {
  layers <- vector("list", length(predictors))
  for (i in seq_along(predictors)) {
    name <- predictors[i]
    if (name %in% names(data)) {
      layers[[i]] <- data[[name]]
    } else if (!is.null(onto_grid) && name %in% names(onto_grid)) {
      layers[[i]] <- terra::resample(onto_grid[[name]], response_template, method = aggregate)
    } else if (is.null(onto_grid)) {
      stop(sprintf(paste0(
        "predictor `%s` named in `formula` is not a layer of `data`.\n",
        "  data layers: %s\n",
        "With an sf/SpatVector `onto`, every predictor must be a layer of `data`; ",
        "the derive-from-onto shortcut needs a gridded `onto`."),
        name, paste(names(data), collapse = ", ")))
    } else {
      stop(sprintf(paste0(
        "predictor `%s` named in `formula` is in neither `data` nor `onto`.\n",
        "  data layers: %s\n  onto layers: %s"),
        name, paste(names(data), collapse = ", "), paste(names(onto_grid), collapse = ", ")))
    }
    names(layers[[i]]) <- name
  }
  terra::rast(layers)
}

# Fit the coarse coefficient grids. The predictors are shared across responses, so
# one engine call returns the intercept, slopes, and per-cell R-squared/residual SD
# for every response, plus a valid-cell count shared across them (the mask is
# complete-case). Returns the shared predictor names, a per-response fit (the
# coefficient SpatRaster, the R-squared grid, and the residual-SD grid, each on the
# coarse grid), and the shared valid-cell-count grid.
fit_windows <- function(parsed, data, onto_grid, radius, aggregate, min_cells, min_variance) {
  response_raster  <- select_layers(data, parsed$response, "data")
  template         <- response_raster[[1L]]
  predictor_raster <- resolve_predictors(parsed$predictors, data, onto_grid,
                                         template, aggregate)

  response_matrices <- stats::setNames(
    lapply(seq_len(terra::nlyr(response_raster)),
           function(i) raster_to_matrix(response_raster[[i]])),
    parsed$response)
  predictor_matrices <- lapply(seq_len(terra::nlyr(predictor_raster)),
                               function(i) raster_to_matrix(predictor_raster[[i]]))
  regression <- window_regression(response_matrices, predictor_matrices,
                                  radius = radius, min_cells = min_cells,
                                  min_variance = min_variance)

  n_valid <- matrix_to_raster(regression$n_valid, template)
  names(n_valid) <- "n.valid"

  responses <- lapply(parsed$response, function(resp) {
    coefficients <- terra::rast(c(
      list(matrix_to_raster(regression$intercept[[resp]], template)),
      lapply(regression$slope[[resp]], matrix_to_raster, template = template)))
    names(coefficients) <- c("(Intercept)", parsed$predictors)
    r_squared <- matrix_to_raster(regression$r_squared[[resp]], template)
    names(r_squared) <- "r.squared"
    residual_sd <- matrix_to_raster(regression$residual_sd[[resp]], template)
    names(residual_sd) <- "residual.sd"
    list(response = resp, coefficients = coefficients,
         r_squared = r_squared, residual_sd = residual_sd)
  })
  list(predictors = parsed$predictors, responses = responses, n_valid = n_valid)
}

# Evaluate one response's fitted relationship on the target. The coarse coefficient
# grids, the R-squared grid, and the residual-SD grid are brought onto the target
# together (resampled to the fine grid, or interpolated at point geometries) and
# combined with the target predictors as fitted = intercept + sum_j slope_j *
# predictor_j. Returns the fitted values, the per-coefficient grids/columns, the
# R-squared, and the residual SD, in the target's representation.
cast_onto <- function(rfit, predictors, target, method) {
  bundle <- bring_onto_target(c(rfit$coefficients, rfit$r_squared, rfit$residual_sd),
                              target, method)
  preds  <- target_predictors(target, predictors)

  intercept <- bundle[["(Intercept)"]]
  slopes <- stats::setNames(lapply(predictors, function(p) bundle[[p]]), predictors)
  predictor_values <- stats::setNames(lapply(predictors, function(p) preds[[p]]),
                                       predictors)
  fitted <- eval_fitted(intercept, slopes, predictor_values)

  coef_cols <- stats::setNames(
    c(list(intercept), lapply(predictors, function(p) slopes[[p]])),
    c("(Intercept)", predictors))
  list(fitted = fitted, coef_cols = coef_cols,
       r_squared   = clamp_unit(bundle[["r.squared"]]),
       residual_sd = clamp_nonneg(bundle[["residual.sd"]]))
}

# Bring the shared valid-cell-count grid onto the target the same way as any other
# coarse diagnostic. It is reported once per topocast() call, not once per response,
# because the mask that produced it is complete-case across every response fit in
# that call.
cast_n_valid <- function(n_valid, target, method) {
  clamp_nonneg(as_value(bring_onto_target(n_valid, target, method)))
}

# Carry each coarse period's anomaly, relative to a coarse baseline, onto the fine
# baseline. Ratio for non-negative variables, additive otherwise. The baseline and
# periods are brought onto the target the same way as the fit. Returns a named list
# of per-period grids/columns.
carry_anomalies <- function(fine_baseline, anomaly, baseline, target, type, method) {
  coarse_baseline <- as_value(bring_onto_target(baseline, target, method))

  n_periods <- terra::nlyr(anomaly)
  periods <- vector("list", n_periods)
  for (period in seq_len(n_periods)) {
    coarse_value <- as_value(bring_onto_target(anomaly[[period]], target, method))
    if (type == "ratio") {
      periods[[period]] <- fine_baseline * safe_ratio(coarse_value, coarse_baseline)
    } else {
      periods[[period]] <- fine_baseline + (coarse_value - coarse_baseline)
    }
  }
  stats::setNames(periods, names(anomaly))
}

# Raster to a wide matrix (row i, column j of the matrix is raster row i, col j).
raster_to_matrix <- function(raster) terra::as.matrix(raster, wide = TRUE)

# Wide matrix back to a raster on a template grid, in row-major cell order.
matrix_to_raster <- function(matrix, template) {
  terra::setValues(template, as.vector(t(matrix)))
}
