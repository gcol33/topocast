#' Moving-window linear regression over raster matrices
#'
#' Fits, for every cell, a linear regression of one or more responses on a shared
#' set of predictors over a square moving window, and returns full-resolution
#' intercept and slope grids together with a per-cell coefficient of determination.
#' The fit uses summed-area tables, so the cost does not grow with the window size.
#' Cells are excluded from a window where any response or predictor is non-finite; a
#' cell is returned as `NA` when its window holds fewer valid cells than the model
#' needs or a predictor has no spread. If every cell of every response and predictor
#' is non-finite, there is nothing to regress and the call errors instead of
#' returning an all-`NA` result.
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
#' @param radius Integer window radius in cells, at least `1`; the window is a
#'   square of side `2 * radius + 1`. A radius of `0` is rejected: its one-cell
#'   window can never hold the `k + 1` valid cells a fit needs.
#' @param min_cells Integer, additional valid cells required in a window beyond
#'   the `k + 1` model terms (`k` predictors plus the intercept). Default `0`.
#' @param min_variance Numeric, the minimum within-window variance a predictor
#'   must have for the cell to be fit. Default `1e-8`.
#'
#' @return A list with `intercept` (a numeric matrix), `slope` (a list of numeric
#'   matrices, one per predictor), `r_squared` (a numeric matrix), `residual_sd` (a
#'   numeric matrix, the residual standard deviation of the local fit; `NA` where
#'   the window has no residual degrees of freedom), and `n_valid` (a numeric
#'   matrix, the count of valid cells the window held), each the same size as `y`.
#'   `n_valid` is shared across responses, since the valid-cell mask is
#'   complete-case across them. When `y` is a list of responses, `intercept`,
#'   `r_squared`, and `residual_sd` are lists of matrices and `slope` is a list of
#'   per-response slope lists, named for the responses.
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
  if (!single) {
    nm <- names(y_list)
    partially_named <- !is.null(nm) && any(nm == "")
    if (partially_named)
      stop(sprintf(paste0(
        "every element of `y` must be named, or none of them; got unnamed element(s) ",
        "at position(s) %s. A result for an unnamed response is not reachable by name."),
        paste(which(nm == ""), collapse = ", ")))
    # A repeated name has the same "unreachable by name" problem as an unnamed
    # element: `res$slope[["name"]]` always resolves to the first match, so the
    # duplicate's result would silently vanish rather than erroring.
    if (!is.null(nm) && anyDuplicated(nm))
      stop(sprintf(paste0(
        "every element of `y` must have a unique name; got a repeated name: %s. ",
        "A result under a repeated name is not reachable by name."),
        paste(unique(nm[duplicated(nm)]), collapse = ", ")))
  }

  radius       <- check_count(radius, "radius")
  min_cells    <- check_count(min_cells, "min_cells")
  min_variance <- check_nonneg(min_variance, "min_variance")

  # A radius of 0 gives a window of exactly one cell, which can never hold the
  # k + 1 valid cells a fit with k predictors needs (k >= 1 is already required
  # above), so every cell of the result would silently be NA. Reject it rather
  # than let that happen with no error or warning.
  if (radius < 1L)
    stop(sprintf(paste0(
      "`radius` must be at least 1; a radius of 0 gives a window of a single cell, ",
      "which can never hold the %d valid cells needed to fit an intercept plus ",
      "%d predictor(s)."), length(x) + 1L, length(x)))

  # A radius at or beyond the grid dimensions already covers every cell in
  # every window; clamping here keeps `row + radius` from overflowing a signed
  # 32-bit int in the C++ engine for a radius near .Machine$integer.max, without
  # changing any result.
  grid_span <- max(nrow(y_list[[1L]]), ncol(y_list[[1L]]))
  radius <- min(radius, grid_span)

  res <- tryCatch(
    window_regression_cpp(y_list, x, radius, min_cells, min_variance),
    error = function(e) {
      # Translate the C++ engine's generic message into one that names the
      # actual R arguments, rather than letting the raw C++ text through.
      if (grepl("no finite cells to regress", conditionMessage(e), fixed = TRUE))
        stop("`y` and `x` share no cell where every response and every predictor ",
             "is finite; there is nothing to regress.", call. = FALSE)
      stop(e)
    })

  if (single)
    return(list(intercept   = res$intercept[[1L]],
                slope       = res$slope[[1L]],
                r_squared   = res$r_squared[[1L]],
                residual_sd = res$residual_sd[[1L]],
                n_valid     = res$n_valid))

  nm <- names(y_list)
  if (is.null(nm)) nm <- paste0("response", seq_along(y_list))
  list(intercept   = stats::setNames(res$intercept, nm),
       slope       = stats::setNames(res$slope, nm),
       r_squared   = stats::setNames(res$r_squared, nm),
       residual_sd = stats::setNames(res$residual_sd, nm),
       n_valid     = res$n_valid)
}

# A single finite, non-negative whole number, coerced to integer. `radius` and
# `min_cells` reach the C++ engine's window bounds and valid-cell threshold
# unchecked otherwise; a negative radius there indexes the summed-area table out of
# bounds and crashes the R session rather than erroring. The upper bound matters as
# much as the lower one: a value above `.Machine$integer.max` still passes a
# "non-negative whole number" test but silently becomes `NA` on `as.integer()`,
# which reaches the C++ engine as `INT_MIN` and crashes the session the same way.
check_count <- function(x, argument) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0 && x == round(x) &&
    x <= .Machine$integer.max
  if (!ok)
    stop(sprintf("`%s` must be a single non-negative whole number no larger than %s, not %s",
                 argument, .Machine$integer.max, paste(deparse(x), collapse = " ")))
  as.integer(x)
}

# A single finite, non-negative number. `min_variance` reaches the C++ engine's
# no-spread guard unchecked otherwise; `NA` or a negative value compares as false
# against any window variance, silently disabling the guard rather than rejecting a
# degenerate predictor as documented.
check_nonneg <- function(x, argument) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0
  if (!ok)
    stop(sprintf("`%s` must be a single non-negative number, not %s",
                 argument, paste(deparse(x), collapse = " ")))
  as.numeric(x)
}
