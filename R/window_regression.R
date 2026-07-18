#' Moving-window linear regression over raster matrices
#'
#' Fits, for every cell, a linear regression of one or more responses on a shared
#' set of predictors over a square moving window, and returns full-resolution
#' intercept and slope grids together with a per-cell coefficient of determination.
#' The fit uses summed-area tables, so the cost does not grow with the window size.
#' Cells are excluded from a window where any response or predictor is non-finite; a
#' cell is returned as `NA` when its window holds fewer valid cells than the model
#' needs or a predictor has no spread.
#'
#' Several responses share the predictors, so the window design is assembled and
#' factored once and solved against every response, and each extra response costs
#' only a back-substitution. The valid-cell mask is therefore complete-case across
#' the responses: a cell enters a window only where every response and predictor is
#' finite.
#'
#' This is the matrix engine behind [topocast()]; it works on plain numeric
#' matrices and does not depend on terra.
#'
#' @param y Numeric matrix, the response on the coarse grid, or a list of numeric
#'   matrices for several responses sharing the predictors and grid.
#' @param x Numeric matrix, or a list of numeric matrices, the predictor(s) on the
#'   same grid as `y`.
#' @param radius Integer window radius in cells; the window is a square of side
#'   `2 * radius + 1`.
#' @param min_cells Integer, additional valid cells required in a window beyond
#'   the `k + 1` model terms (`k` predictors plus the intercept). Default `0`.
#' @param min_variance Numeric, the minimum within-window variance a predictor
#'   must have for the cell to be fit. Default `1e-8`.
#'
#' @return A list with `intercept` (a numeric matrix), `slope` (a list of numeric
#'   matrices, one per predictor), and `r_squared` (a numeric matrix), each the
#'   same size as `y`. When `y` is a list of responses, `intercept` and `r_squared`
#'   are lists of matrices and `slope` is a list of per-response slope lists, named
#'   for the responses.
#'
#' @seealso [topocast()] for the terra workflow.
#'
#' @examples
#' set.seed(1)
#' elevation <- matrix(runif(100, 0, 1000), 10, 10)
#' climate <- 50 + 0.01 * elevation
#' fit <- window_regression(climate, elevation, radius = 3)
#' str(fit)
#'
#' # several responses sharing the predictor: the design is factored once
#' rain <- 800 - 0.1 * elevation
#' temp <- 15 - 0.006 * elevation
#' both <- window_regression(list(rain = rain, temp = temp), elevation, radius = 3)
#' names(both$slope)
#'
#' @export
window_regression <- function(y, x, radius, min_cells = 0L, min_variance = 1e-8) {
  if (is.matrix(x)) x <- list(x)
  if (!is.list(x) || length(x) < 1L)
    stop("`x` must be a matrix or a list of matrices")
  if (!all(vapply(x, is.matrix, logical(1))))
    stop("every element of `x` must be a matrix")

  single <- is.matrix(y)
  y_list <- if (single) list(y) else y
  if (!is.list(y_list) || length(y_list) < 1L)
    stop("`y` must be a matrix or a list of matrices")
  if (!all(vapply(y_list, is.matrix, logical(1))))
    stop("every element of `y` must be a matrix")

  radius     <- check_count(radius, "radius")
  min_cells  <- check_count(min_cells, "min_cells")

  res <- window_regression_cpp(y_list, x, radius, min_cells, as.numeric(min_variance))

  if (single)
    return(list(intercept = res$intercept[[1L]],
                slope     = res$slope[[1L]],
                r_squared = res$r_squared[[1L]]))

  nm <- names(y_list)
  if (is.null(nm)) nm <- paste0("response", seq_along(y_list))
  list(intercept = stats::setNames(res$intercept, nm),
       slope     = stats::setNames(res$slope, nm),
       r_squared = stats::setNames(res$r_squared, nm))
}

# A single finite, non-negative whole number, coerced to integer. `radius` and
# `min_cells` reach the C++ engine's window bounds and valid-cell threshold
# unchecked otherwise; a negative radius there indexes the summed-area table out of
# bounds and crashes the R session rather than erroring.
check_count <- function(x, argument) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0 && x == round(x)
  if (!ok)
    stop(sprintf("`%s` must be a single non-negative whole number, not %s",
                 argument, paste(deparse(x), collapse = " ")))
  as.integer(x)
}
