# Brute-force reference: an independent per-cell OLS over the window, used only as
# a test oracle for the summed-area-table engine.
window_regression_ref <- function(y, x, radius, min_cells = 0L, min_variance = 1e-8) {
  if (is.matrix(x)) x <- list(x)
  rows <- nrow(y); cols <- ncol(y); k <- length(x)
  intercept <- matrix(NA_real_, rows, cols)
  r_squared <- matrix(NA_real_, rows, cols)
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
    fit <- stats::lm.fit(cbind(1, x_valid), y_valid)
    intercept[row, col] <- fit$coefficients[1]
    for (j in seq_len(k)) slope[[j]][row, col] <- fit$coefficients[1 + j]
    sst <- sum((y_valid - mean(y_valid))^2)
    if (sst > 0) r_squared[row, col] <- max(0, min(1, 1 - sum(fit$residuals^2) / sst))
  }
  list(intercept = intercept, slope = slope, r_squared = r_squared)
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

test_that("exact linear data has r_squared 1, noisy data matches the oracle", {
  set.seed(21)
  elevation <- matrix(runif(144, 0, 2000), 12, 12)
  exact <- 30 - 0.006 * elevation
  fit_exact <- window_regression(exact, elevation, radius = 4)
  finite <- is.finite(fit_exact$r_squared)
  expect_true(all(fit_exact$r_squared[finite] > 1 - 1e-8))

  x1 <- matrix(runif(150, 0, 2500), 10, 15)
  x2 <- matrix(runif(150, 0, 50), 10, 15)
  y <- 5 - 0.004 * x1 + 0.2 * x2 + matrix(rnorm(150, 0, 0.3), 10, 15)
  fit <- window_regression(y, list(x1, x2), radius = 4)
  ref <- window_regression_ref(y, list(x1, x2), radius = 4)
  expect_equal(fit$r_squared, ref$r_squared, tolerance = 1e-6)
})

test_that("several responses match per-response single fits and are named", {
  set.seed(22)
  elevation <- matrix(runif(120, 0, 3000), 10, 12)
  rain <- 800 - 0.1 * elevation + matrix(rnorm(120, 0, 5), 10, 12)
  temp <- 15 - 0.006 * elevation + matrix(rnorm(120, 0, 0.3), 10, 12)

  multi <- window_regression(list(rain = rain, temp = temp), elevation, radius = 3)
  expect_equal(names(multi$slope), c("rain", "temp"))
  expect_equal(names(multi$intercept), c("rain", "temp"))

  one_rain <- window_regression(rain, elevation, radius = 3)
  one_temp <- window_regression(temp, elevation, radius = 3)
  expect_equal(multi$slope[["rain"]][[1]], one_rain$slope[[1]])
  expect_equal(multi$slope[["temp"]][[1]], one_temp$slope[[1]])
  expect_equal(multi$r_squared[["temp"]], one_temp$r_squared)
})

test_that("the complete-case mask is shared across responses", {
  set.seed(23)
  elevation <- matrix(runif(120, 0, 3000), 10, 12)
  rain <- 800 - 0.1 * elevation + matrix(rnorm(120, 0, 5), 10, 12)
  temp <- 15 - 0.006 * elevation + matrix(rnorm(120, 0, 0.3), 10, 12)
  temp[4, 5] <- NA   # missing in one response only

  multi <- window_regression(list(rain, temp), elevation, radius = 2)
  # the shared mask drops cell (4,5) from rain's fit too: rain alone keeps it.
  rain_masked <- rain; rain_masked[4, 5] <- NA
  ref <- window_regression_ref(rain_masked, elevation, radius = 2)
  expect_equal(multi$slope[[1]][[1]], ref$slope[[1]], tolerance = 1e-6)
})

test_that("radius and min_cells are validated as non-negative whole numbers (issue #4)", {
  set.seed(24)
  elevation <- matrix(runif(100, 0, 2000), 10, 10)
  climate <- 30 - 0.006 * elevation

  # a negative radius previously indexed the summed-area table out of bounds and
  # crashed the R session; it must now be a clean error.
  expect_error(window_regression(climate, elevation, radius = -3), "non-negative whole number")
  expect_error(window_regression(climate, elevation, radius = 2.9), "non-negative whole number")
  expect_error(window_regression(climate, elevation, radius = NA), "non-negative whole number")
  expect_error(window_regression(climate, elevation, radius = c(3, 4)), "non-negative whole number")
  expect_error(window_regression(climate, elevation, radius = 3, min_cells = -5),
               "non-negative whole number")

  # valid whole-number radii, integer or double, still work and agree
  expect_equal(window_regression(climate, elevation, radius = 3),
               window_regression(climate, elevation, radius = 3L))
})

test_that("a radius/min_cells above .Machine$integer.max errors instead of crashing (issue #9)", {
  set.seed(24)
  elevation <- matrix(runif(100, 0, 2000), 10, 10)
  climate <- 30 - 0.006 * elevation

  # a radius/min_cells too large for a 32-bit int previously overflowed to NA on
  # as.integer() and crashed the R session; it must now be a clean error.
  expect_error(window_regression(climate, elevation, radius = 3e9), "non-negative whole number")
  expect_error(window_regression(climate, elevation, radius = 3, min_cells = 3e9),
               "non-negative whole number")

  # one past the boundary is still rejected; .Machine$integer.max itself is a
  # valid (if impractical) radius and is not exercised here.
  expect_error(window_regression(climate, elevation, radius = .Machine$integer.max + 1),
               "non-negative whole number")
})

test_that("min_variance is validated as a non-negative number (issue #10)", {
  set.seed(26)
  flat <- matrix(1500, 8, 8)
  climate <- matrix(rnorm(64, 10, 1), 8, 8)

  # NA/negative min_variance previously compared as false against every window
  # variance, silently disabling the no-spread guard instead of rejecting it.
  expect_error(window_regression(climate, flat, radius = 3, min_variance = NA),
               "non-negative number")
  expect_error(window_regression(climate, flat, radius = 3, min_variance = -1),
               "non-negative number")
  expect_error(window_regression(climate, flat, radius = 3, min_variance = c(1e-8, 2e-8)),
               "non-negative number")

  # the guard still holds at the default: a flat predictor yields NA, not a
  # plausible-looking coefficient from a singular design.
  fit <- window_regression(climate, flat, radius = 3)
  expect_true(all(is.na(fit$slope[[1]])))
  expect_true(all(is.na(fit$intercept)))
})

test_that("a partially-named response list is rejected rather than losing a result (issue #6)", {
  set.seed(25)
  elevation <- matrix(runif(100, 0, 2000), 10, 10)
  rain <- 800 - 0.1 * elevation
  temp <- 15 - 0.006 * elevation

  expect_error(window_regression(list(rain = rain, temp), elevation, radius = 3),
               "must be named, or none of them")

  # fully unnamed and fully named lists both still work
  unnamed <- window_regression(list(rain, temp), elevation, radius = 3)
  expect_equal(names(unnamed$slope), c("response1", "response2"))
  named <- window_regression(list(rain = rain, temp = temp), elevation, radius = 3)
  expect_equal(names(named$slope), c("rain", "temp"))
})
