#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>

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

// Moving-window linear regression of a response on k predictors. For every cell
// the regression is fit over the square window of the given radius (in cells) and
// the fitted intercept and slopes are returned as full-resolution grids. Response
// and predictors are centred on their global means before accumulation to keep the
// summed-area tables well conditioned; the intercept is returned on the raw scale.
//
// [[Rcpp::export]]
List window_regression_cpp(const arma::mat& y, const List& Xlist, int radius,
                           int min_cells, double min_variance) {
  const arma::uword rows = y.n_rows, cols = y.n_cols;
  const int k = Xlist.size();
  if (k < 1) stop("at least one predictor is required");

  std::vector<arma::mat> X;
  X.reserve(k);
  for (int j = 0; j < k; ++j) {
    arma::mat xj = as<arma::mat>(Xlist[j]);
    if (xj.n_rows != rows || xj.n_cols != cols)
      stop("a predictor grid does not match the response dimensions");
    X.push_back(xj);
  }

  // A cell contributes only where the response and every predictor are finite.
  arma::mat mask(rows, cols, arma::fill::ones);
  for (arma::uword i = 0; i < rows * cols; ++i) {
    bool ok = std::isfinite(y(i));
    for (int j = 0; j < k && ok; ++j) ok = std::isfinite(X[j](i));
    if (!ok) mask(i) = 0.0;
  }

  const double n_all = arma::accu(mask);
  if (n_all < 1) stop("no finite cells to regress");

  // Global means over contributing cells, used to centre before accumulation.
  double mean_y = 0.0;
  std::vector<double> mean_x(k, 0.0);
  for (arma::uword i = 0; i < rows * cols; ++i) {
    if (mask(i) > 0) {
      mean_y += y(i);
      for (int j = 0; j < k; ++j) mean_x[j] += X[j](i);
    }
  }
  mean_y /= n_all;
  for (int j = 0; j < k; ++j) mean_x[j] /= n_all;

  // Centred, masked copies (invalid cells become 0 and so drop out of the sums).
  arma::mat y_centred(rows, cols, arma::fill::zeros);
  std::vector<arma::mat> X_centred(k, arma::mat(rows, cols, arma::fill::zeros));
  for (arma::uword i = 0; i < rows * cols; ++i) {
    if (mask(i) > 0) {
      y_centred(i) = y(i) - mean_y;
      for (int j = 0; j < k; ++j) X_centred[j](i) = X[j](i) - mean_x[j];
    }
  }

  // Summed-area tables for the sufficient statistics of the normal equations: the
  // count, the right-hand side X'y, and the cross-products X'X.
  arma::mat ii_count = integral_image(mask);
  arma::mat ii_y = integral_image(y_centred);
  std::vector<arma::mat> ii_x(k), ii_xy(k);
  std::vector<std::vector<arma::mat>> ii_xx(k, std::vector<arma::mat>(k));
  for (int j = 0; j < k; ++j) {
    ii_x[j]  = integral_image(X_centred[j]);
    ii_xy[j] = integral_image(X_centred[j] % y_centred);
    for (int l = j; l < k; ++l) ii_xx[j][l] = integral_image(X_centred[j] % X_centred[l]);
  }

  arma::mat intercept(rows, cols);
  intercept.fill(arma::datum::nan);
  arma::mat nan_grid(rows, cols);
  nan_grid.fill(arma::datum::nan);
  std::vector<arma::mat> slope(k, nan_grid);

  const int need = k + 1 + (min_cells > 0 ? min_cells : 0);

  for (arma::uword row = 0; row < rows; ++row) {
    const int r0 = std::max(0, (int)row - radius);
    const int r1 = std::min((int)rows - 1, (int)row + radius);
    for (arma::uword col = 0; col < cols; ++col) {
      const int c0 = std::max(0, (int)col - radius);
      const int c1 = std::min((int)cols - 1, (int)col + radius);

      const double n = window_sum(ii_count, r0, c0, r1, c1);
      if (n < need) continue;

      arma::mat A(k + 1, k + 1, arma::fill::zeros);
      arma::vec b(k + 1, arma::fill::zeros);
      A(0, 0) = n;
      b(0) = window_sum(ii_y, r0, c0, r1, c1);

      bool degenerate = false;
      for (int j = 0; j < k; ++j) {
        const double sum_x = window_sum(ii_x[j], r0, c0, r1, c1);
        A(0, j + 1) = sum_x;
        A(j + 1, 0) = sum_x;
        b(j + 1) = window_sum(ii_xy[j], r0, c0, r1, c1);
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

      arma::vec beta;
      if (!arma::solve(beta, A, b)) continue;

      // Map the centred intercept back onto the raw response/predictor scale.
      double a = beta(0) + mean_y;
      for (int j = 0; j < k; ++j) {
        slope[j](row, col) = beta(j + 1);
        a -= beta(j + 1) * mean_x[j];
      }
      intercept(row, col) = a;
    }
  }

  List slope_list(k);
  for (int j = 0; j < k; ++j) slope_list[j] = slope[j];
  return List::create(_["intercept"] = intercept, _["slope"] = slope_list);
}
