// Armadillo warnings are silenced so the per-cell solve never writes to an R
// stream from inside the parallel region below, where concurrent stream access
// would be unsafe. The solve still reports failure through its return value.
#define ARMA_WARN_LEVEL 0

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Summed-area table with a zero-padded first row and column, so the sum over any
// rectangular window is four lookups regardless of the window size. Invalid cells
// must already carry 0 (they then contribute nothing to the sums; their count is
// tracked by the separate mask table).
static arma::mat integral_image(const arma::mat& data) {
  const arma::uword rows = data.n_rows, cols = data.n_cols;
  arma::mat ii(rows + 1, cols + 1, arma::fill::zeros);
  for (arma::uword c = 0; c < cols; ++c) {
    for (arma::uword r = 0; r < rows; ++r) {
      ii(r + 1, c + 1) = data(r, c) + ii(r, c + 1) + ii(r + 1, c) - ii(r, c);
    }
  }
  return ii;
}

// Inclusive sum over rows [r0, r1] and columns [c0, c1] from a padded table.
static inline double window_sum(const arma::mat& ii,
                                int r0, int c0, int r1, int c1) {
  return ii(r1 + 1, c1 + 1) - ii(r0, c1 + 1) - ii(r1 + 1, c0) + ii(r0, c0);
}

// Moving-window linear regression of one or more responses on k shared predictors.
// For every cell the regression is fit over the square window of the given radius
// (in cells); the fitted intercept and slopes are returned as full-resolution grids.
//
// All responses share the predictors, so the window design matrix A = Z'Z depends
// only on the predictors and the valid-cell mask. It is assembled and factored once
// per window and solved against every response at once (the right-hand sides form
// the columns of B), so each extra response adds only a back-substitution. The mask
// is therefore complete-case across the response stack: a cell contributes to a
// window only where every response and every predictor is finite.
//
// Response and predictors are centred on their global means before accumulation to
// keep the summed-area tables well conditioned; the intercept is returned on the raw
// scale. A per-window coefficient of determination and residual standard deviation
// are returned for each response, and a valid-cell count (shared across responses,
// since the mask is complete-case), all from the same sufficient statistics.
//
// [[Rcpp::export]]
List window_regression_cpp(const List& Ylist, const List& Xlist, int radius,
                           int min_cells, double min_variance) {
  const int k = Xlist.size();
  if (k < 1) stop("at least one predictor is required");
  const int R = Ylist.size();
  if (R < 1) stop("at least one response is required");

  std::vector<arma::mat> Y;
  Y.reserve(R);
  arma::uword rows = 0, cols = 0;
  for (int r = 0; r < R; ++r) {
    arma::mat yr = as<arma::mat>(Ylist[r]);
    if (r == 0) { rows = yr.n_rows; cols = yr.n_cols; }
    else if (yr.n_rows != rows || yr.n_cols != cols)
      stop("all response grids must share the same dimensions");
    Y.push_back(yr);
  }

  std::vector<arma::mat> X;
  X.reserve(k);
  for (int j = 0; j < k; ++j) {
    arma::mat xj = as<arma::mat>(Xlist[j]);
    if (xj.n_rows != rows || xj.n_cols != cols)
      stop("a predictor grid does not match the response dimensions");
    X.push_back(xj);
  }

  // A cell contributes only where every response and every predictor is finite.
  arma::mat mask(rows, cols, arma::fill::ones);
  for (arma::uword i = 0; i < rows * cols; ++i) {
    bool ok = true;
    for (int r = 0; r < R && ok; ++r) ok = std::isfinite(Y[r](i));
    for (int j = 0; j < k && ok; ++j) ok = std::isfinite(X[j](i));
    if (!ok) mask(i) = 0.0;
  }

  const double n_all = arma::accu(mask);
  if (n_all < 1) stop("no finite cells to regress");

  // Global means over contributing cells, used to centre before accumulation.
  std::vector<double> mean_y(R, 0.0), mean_x(k, 0.0);
  for (arma::uword i = 0; i < rows * cols; ++i) {
    if (mask(i) > 0) {
      for (int r = 0; r < R; ++r) mean_y[r] += Y[r](i);
      for (int j = 0; j < k; ++j) mean_x[j] += X[j](i);
    }
  }
  for (int r = 0; r < R; ++r) mean_y[r] /= n_all;
  for (int j = 0; j < k; ++j) mean_x[j] /= n_all;

  // Centred, masked copies (invalid cells become 0 and so drop out of the sums).
  std::vector<arma::mat> Y_c(R, arma::mat(rows, cols, arma::fill::zeros));
  std::vector<arma::mat> X_c(k, arma::mat(rows, cols, arma::fill::zeros));
  for (arma::uword i = 0; i < rows * cols; ++i) {
    if (mask(i) > 0) {
      for (int r = 0; r < R; ++r) Y_c[r](i) = Y[r](i) - mean_y[r];
      for (int j = 0; j < k; ++j) X_c[j](i) = X[j](i) - mean_x[j];
    }
  }

  // Summed-area tables for the sufficient statistics of the normal equations: the
  // count, the predictor sums and cross-products X'X (shared across responses), and
  // per response the right-hand side X'y plus y'y for the coefficient of determination.
  arma::mat ii_count = integral_image(mask);
  std::vector<arma::mat> ii_y(R), ii_yy(R);
  for (int r = 0; r < R; ++r) {
    ii_y[r]  = integral_image(Y_c[r]);
    ii_yy[r] = integral_image(Y_c[r] % Y_c[r]);
  }
  std::vector<arma::mat> ii_x(k);
  std::vector<std::vector<arma::mat>> ii_xx(k, std::vector<arma::mat>(k));
  std::vector<std::vector<arma::mat>> ii_xy(k, std::vector<arma::mat>(R));
  for (int j = 0; j < k; ++j) {
    ii_x[j] = integral_image(X_c[j]);
    for (int r = 0; r < R; ++r) ii_xy[j][r] = integral_image(X_c[j] % Y_c[r]);
    for (int l = j; l < k; ++l) ii_xx[j][l] = integral_image(X_c[j] % X_c[l]);
  }

  arma::mat nan_grid(rows, cols);
  nan_grid.fill(arma::datum::nan);
  std::vector<arma::mat> intercept(R, nan_grid);
  std::vector<std::vector<arma::mat>> slope(R, std::vector<arma::mat>(k, nan_grid));
  std::vector<arma::mat> r_squared(R, nan_grid);
  std::vector<arma::mat> residual_sd(R, nan_grid);
  // The valid-cell mask is shared across every response in this call (see the
  // complete-case masking above), so one grid of window counts covers all of them.
  arma::mat n_valid = nan_grid;

  // Computed in double, not int: k + 1 + min_cells can exceed INT_MAX when
  // min_cells is large (it is validated only against .Machine$integer.max on
  // the R side, not against k + 1 + min_cells), and n below is already a double.
  const double need = (double)k + 1.0 + (double)min_cells;
  const int nrows = (int)rows, ncols = (int)cols;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int row = 0; row < nrows; ++row) {
    const int r0 = std::max(0, row - radius);
    const int r1 = std::min(nrows - 1, row + radius);
    for (int col = 0; col < ncols; ++col) {
      const int c0 = std::max(0, col - radius);
      const int c1 = std::min(ncols - 1, col + radius);

      const double n = window_sum(ii_count, r0, c0, r1, c1);
      if (n < need) continue;

      // Recorded as soon as the window holds enough valid cells to attempt a
      // fit, independent of whether the fit itself succeeds below: n_valid
      // reports how much data the window held, not whether that data produced
      // a coefficient. r_squared/residual_sd/intercept/slope stay NA in the
      // degenerate/singular cases that follow.
      n_valid(row, col) = n;

      // Predictor side of the normal equations, shared by every response.
      arma::mat A(k + 1, k + 1, arma::fill::zeros);
      A(0, 0) = n;
      bool degenerate = false;
      for (int j = 0; j < k; ++j) {
        const double sum_x = window_sum(ii_x[j], r0, c0, r1, c1);
        A(0, j + 1) = sum_x;
        A(j + 1, 0) = sum_x;
        const double sum_xx = window_sum(ii_xx[j][j], r0, c0, r1, c1);
        A(j + 1, j + 1) = sum_xx;
        // Reject windows where a predictor has effectively no spread.
        const double variance = sum_xx / n - (sum_x / n) * (sum_x / n);
        if (variance < min_variance) { degenerate = true; break; }
        for (int l = j + 1; l < k; ++l) {
          const double sum_xl = window_sum(ii_xx[j][l], r0, c0, r1, c1);
          A(j + 1, l + 1) = sum_xl;
          A(l + 1, j + 1) = sum_xl;
        }
      }
      if (degenerate) continue;

      // Response side: one right-hand-side column per response.
      arma::mat B(k + 1, R, arma::fill::zeros);
      for (int r = 0; r < R; ++r) {
        B(0, r) = window_sum(ii_y[r], r0, c0, r1, c1);
        for (int j = 0; j < k; ++j)
          B(j + 1, r) = window_sum(ii_xy[j][r], r0, c0, r1, c1);
      }

      // solve_opts::no_approx is required here: Armadillo's default solve()
      // silently falls back to a least-squares (pseudo-inverse) solution on a
      // detected-singular system instead of failing, which would return a
      // coefficient for an exactly rank-deficient window (e.g. one predictor
      // an exact linear function of another) rather than the documented NA.
      arma::mat Beta;
      if (!arma::solve(Beta, A, B, arma::solve_opts::no_approx)) continue;

      for (int r = 0; r < R; ++r) {
        // Map the centred intercept back onto the raw response/predictor scale.
        double a = Beta(0, r) + mean_y[r];
        for (int j = 0; j < k; ++j) {
          slope[r][j](row, col) = Beta(j + 1, r);
          a -= Beta(j + 1, r) * mean_x[j];
        }
        intercept[r](row, col) = a;

        // R^2 and residual SD from the same sufficient statistics:
        // SSE = y'y - beta'(X'y), SST = y'y - (sum y)^2 / n, both about the window
        // mean of the response.
        const double syy = window_sum(ii_yy[r], r0, c0, r1, c1);
        const double sy  = B(0, r);
        const double sst = syy - sy * sy / n;
        const double sse = syy - arma::dot(Beta.col(r), B.col(r));
        if (sst > 0) {
          double rsq = 1.0 - sse / sst;
          if (rsq < 0.0) rsq = 0.0;
          if (rsq > 1.0) rsq = 1.0;
          r_squared[r](row, col) = rsq;
        }

        // Residual SD needs residual degrees of freedom (n minus the k+1 fitted
        // terms) to be positive; a window with n == need (min_cells == 0) has none.
        const double dof = n - (k + 1);
        if (dof > 0) residual_sd[r](row, col) = std::sqrt(std::max(0.0, sse) / dof);
      }
    }
  }

  List intercept_out(R), slope_out(R), r2_out(R), resid_sd_out(R);
  for (int r = 0; r < R; ++r) {
    intercept_out[r] = intercept[r];
    List slope_r(k);
    for (int j = 0; j < k; ++j) slope_r[j] = slope[r][j];
    slope_out[r] = slope_r;
    r2_out[r] = r_squared[r];
    resid_sd_out[r] = residual_sd[r];
  }
  return List::create(_["intercept"] = intercept_out,
                      _["slope"] = slope_out,
                      _["r_squared"] = r2_out,
                      _["residual_sd"] = resid_sd_out,
                      _["n_valid"] = n_valid);
}
