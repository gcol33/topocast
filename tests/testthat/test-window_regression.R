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

test_that("a duplicate-named response list is rejected rather than hiding a result (issue #15)", {
  set.seed(27)
  elevation <- matrix(runif(100, 0, 2000), 10, 10)
  jan <- 30 - 0.006 * elevation
  jul <- 10 + 0.01 * elevation   # a genuinely different response, same name by mistake

  # previously: no error, and jul's fit was computed but unreachable under
  # `$slope[["prec"]]`, which always resolved to jan's (the first match).
  expect_error(window_regression(list(prec = jan, prec = jul), elevation, radius = 2),
               "unique name")
})

test_that("min_cells near .Machine$integer.max no longer overflows the valid-cell guard (issue #13)", {
  set.seed(28)
  elevation <- matrix(runif(25, 0, 2000), 5, 5)
  climate <- 30 - 0.006 * elevation

  # previously, k + 1 + min_cells overflowed a 32-bit int and wrapped negative,
  # so every window passed the valid-cell check regardless of how few cells it
  # actually held; a 5x5 grid can never hold this many valid cells anywhere.
  fit <- window_regression(climate, elevation, radius = 1,
                           min_cells = .Machine$integer.max - 1)
  expect_true(all(is.na(fit$intercept)))
  expect_true(all(is.na(fit$slope[[1]])))
})

test_that("radius near .Machine$integer.max is clamped to the grid span instead of crashing (issue #14)", {
  set.seed(29)
  elevation <- matrix(runif(15, 0, 2000), 3, 5)
  climate <- 30 - 0.006 * elevation

  # previously, row + radius overflowed a signed 32-bit int for row >= 1 and
  # crashed the R session; any radius at or beyond the grid span covers every
  # cell in every window, so clamping it changes nothing about the result.
  huge <- window_regression(climate, elevation, radius = .Machine$integer.max - 1L)
  full <- window_regression(climate, elevation,
                            radius = max(nrow(elevation), ncol(elevation)))
  expect_equal(huge, full)
})

test_that("n_valid reports the raw valid-cell count even when the window is degenerate (issue #16)", {
  flat <- matrix(1500, 8, 8)
  climate <- matrix(rnorm(64, 10, 1), 8, 8)
  radius <- 3
  fit <- window_regression(climate, flat, radius = radius)

  # the no-spread guard still returns NA coefficients everywhere
  expect_true(all(is.na(fit$slope[[1]])))
  expect_true(all(is.na(fit$intercept)))

  # but n_valid, unlike the coefficients, reports the window's raw valid-cell
  # count rather than NA: every window here holds real data, just a flat
  # predictor, and n_valid should say so rather than looking identical to a
  # window that never had enough cells to begin with.
  expect_false(anyNA(fit$n_valid))
  expect_equal(fit$n_valid[4, 4], (2 * radius + 1)^2)
})

test_that("min_cells > 0 excludes windows that min_cells = 0 would fit (issue #17)", {
  elevation <- matrix(NA_real_, 5, 5)
  climate   <- matrix(NA_real_, 5, 5)
  # exactly two valid, distinct-elevation cells fall in the radius-1 window
  # centred at (3, 3); every other cell in this grid is NA.
  elevation[3, 3] <- 1000; elevation[3, 4] <- 1200
  climate[3, 3]   <- 20;   climate[3, 4]   <- 18

  fit0 <- window_regression(climate, elevation, radius = 1, min_cells = 0)
  expect_false(is.na(fit0$intercept[3, 3]))

  fit1 <- window_regression(climate, elevation, radius = 1, min_cells = 1)
  expect_true(is.na(fit1$intercept[3, 3]))
})

test_that("a rank-deficient multi-predictor window returns NA rather than an unstable fit (issue #17)", {
  # x2 is an exact linear function of x1 (and the intercept) everywhere, so the
  # window design matrix is singular even though x1 and x2 each individually
  # have real within-window variance and so do not trip the no-spread guard.
  x1 <- matrix(1:25, 5, 5)
  x2 <- 3 + 2 * x1
  y  <- matrix(rnorm(25), 5, 5)

  fit <- window_regression(y, list(x1, x2), radius = 2)
  expect_true(all(is.na(fit$intercept)))
  expect_true(all(is.na(fit$slope[[1]])))
  expect_true(all(is.na(fit$slope[[2]])))
})

test_that("a flat window on a large grid still returns NA (issue #21)", {
  # A summed-area-table entry near the far side of a large grid holds a sum over
  # most of the grid's cells, so its magnitude grows with the grid size; a plain
  # double-precision table loses precision proportional to that magnitude, not to
  # the window's own size, and can turn a genuinely zero within-window variance
  # into a small positive number that wrongly clears min_variance. 1200x1200 cells
  # (1.44e6) is comfortably past the point (~1e6) where this was observed to
  # misfire; real DEM tiles (SRTM 1-arcsec: 3601x3601) are larger still.
  set.seed(30)
  n <- 1200
  elevation <- matrix(runif(n * n, 0, 8000), n, n)
  elevation[(n - 6):(n - 2), (n - 6):(n - 2)] <- 4000  # exactly flat 5x5 block
  climate <- matrix(rnorm(n * n, 10, 2), n, n)

  fit <- window_regression(climate, elevation, radius = 2)
  center <- n - 4
  expect_true(is.na(fit$intercept[center, center]))
  expect_true(is.na(fit$slope[[1]][center, center]))
})

test_that("radius = 0 is rejected rather than silently returning an all-NA result (issue #25)", {
  set.seed(31)
  elevation <- matrix(runif(100, 0, 2000), 10, 10)
  climate <- 30 - 0.006 * elevation
  expect_error(window_regression(climate, elevation, radius = 0), "at least 1")
})

test_that("an entirely non-finite y/x errors clearly instead of returning all-NA (issue #30/#34)", {
  y <- matrix(NA_real_, 5, 5)
  x <- matrix(NA_real_, 5, 5)
  expect_error(window_regression(y, x, radius = 2), "nothing to regress")
})

test_that("the engine matches the brute-force oracle with three predictors (issue #30)", {
  set.seed(32)
  x1 <- matrix(runif(120, 0, 2000), 10, 12)
  x2 <- matrix(runif(120, 0, 50), 10, 12)
  x3 <- matrix(runif(120, 0, 10), 10, 12)
  y <- 5 - 0.004 * x1 + 0.2 * x2 + 0.5 * x3 + matrix(rnorm(120, 0, 0.3), 10, 12)

  fit <- window_regression(y, list(x1, x2, x3), radius = 4)
  ref <- window_regression_ref(y, list(x1, x2, x3), radius = 4)
  expect_equal(fit$intercept, ref$intercept, tolerance = 1e-5)
  expect_equal(fit$slope[[1]], ref$slope[[1]], tolerance = 1e-5)
  expect_equal(fit$slope[[2]], ref$slope[[2]], tolerance = 1e-5)
  expect_equal(fit$slope[[3]], ref$slope[[3]], tolerance = 1e-5)
})
