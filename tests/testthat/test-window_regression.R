# Brute-force reference: an independent per-cell OLS over the window, used only as
# a test oracle for the summed-area-table engine.
window_regression_ref <- function(y, x, radius, min_cells = 0L, min_variance = 1e-8) {
  if (is.matrix(x)) x <- list(x)
  rows <- nrow(y); cols <- ncol(y); k <- length(x)
  intercept <- matrix(NA_real_, rows, cols)
  slope <- replicate(k, matrix(NA_real_, rows, cols), simplify = FALSE)
  for (row in seq_len(rows)) for (col in seq_len(cols)) {
    row_range <- max(1, row - radius):min(rows, row + radius)
    col_range <- max(1, col - radius):min(cols, col + radius)
    y_window <- as.vector(y[row_range, col_range, drop = FALSE])
    x_window <- matrix(vapply(x, function(m) as.vector(m[row_range, col_range, drop = FALSE]),
                              numeric(length(y_window))), ncol = k)
    valid <- is.finite(y_window) & apply(matrix(is.finite(x_window), ncol = k), 1, all)
    if (sum(valid) < k + 1 + min_cells) next
    y_valid <- y_window[valid]
    x_valid <- x_window[valid, , drop = FALSE]
    predictor_variance <- apply(x_valid, 2, function(v) mean(v^2) - mean(v)^2)
    if (any(predictor_variance < min_variance)) next
    coef <- stats::lm.fit(cbind(1, x_valid), y_valid)$coefficients
    intercept[row, col] <- coef[1]
    for (j in seq_len(k)) slope[[j]][row, col] <- coef[1 + j]
  }
  list(intercept = intercept, slope = slope)
}

test_that("exact linear data is recovered", {
  set.seed(11)
  elevation <- matrix(runif(144, 0, 2000), 12, 12)
  climate <- 30 - 0.006 * elevation
  fit <- window_regression(climate, elevation, radius = 4)
  expect_equal(fit$slope[[1]], matrix(-0.006, 12, 12), tolerance = 1e-8)
  expect_equal(fit$intercept, matrix(30, 12, 12), tolerance = 1e-6)
})

test_that("the engine matches the brute-force oracle (univariate)", {
  set.seed(12)
  elevation <- matrix(runif(120, 0, 3000), 10, 12)
  climate <- 12 - 0.005 * elevation + matrix(rnorm(120, 0, 0.5), 10, 12)
  fit <- window_regression(climate, elevation, radius = 3)
  ref <- window_regression_ref(climate, elevation, radius = 3)
  expect_equal(fit$intercept, ref$intercept, tolerance = 1e-6)
  expect_equal(fit$slope[[1]], ref$slope[[1]], tolerance = 1e-6)
})

test_that("the engine matches the brute-force oracle (multi-predictor)", {
  set.seed(13)
  x1 <- matrix(runif(150, 0, 2500), 10, 15)
  x2 <- matrix(runif(150, 0, 50), 10, 15)
  y <- 5 - 0.004 * x1 + 0.2 * x2 + matrix(rnorm(150, 0, 0.3), 10, 15)
  fit <- window_regression(y, list(x1, x2), radius = 4)
  ref <- window_regression_ref(y, list(x1, x2), radius = 4)
  expect_equal(fit$intercept, ref$intercept, tolerance = 1e-5)
  expect_equal(fit$slope[[1]], ref$slope[[1]], tolerance = 1e-5)
  expect_equal(fit$slope[[2]], ref$slope[[2]], tolerance = 1e-5)
})

test_that("non-finite cells are excluded, matching the oracle", {
  set.seed(14)
  elevation <- matrix(runif(120, 0, 3000), 10, 12)
  climate <- 12 - 0.005 * elevation + matrix(rnorm(120, 0, 0.5), 10, 12)
  climate[3, 4] <- NA
  elevation[7, 9] <- NA
  fit <- window_regression(climate, elevation, radius = 3)
  ref <- window_regression_ref(climate, elevation, radius = 3)
  expect_equal(fit$intercept, ref$intercept, tolerance = 1e-6)
  expect_equal(fit$slope[[1]], ref$slope[[1]], tolerance = 1e-6)
})

test_that("a predictor with no spread yields NA", {
  flat <- matrix(1500, 8, 8)
  climate <- matrix(rnorm(64, 10, 1), 8, 8)
  fit <- window_regression(climate, flat, radius = 3)
  expect_true(all(is.na(fit$slope[[1]])))
  expect_true(all(is.na(fit$intercept)))
})
