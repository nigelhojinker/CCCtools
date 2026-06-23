// fast_group_summary.cpp
#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

double quantile_type7(std::vector<double> x, double p) {
  int n = x.size();

  if (n == 0) return NA_REAL;
  if (n == 1) return x[0];

  std::sort(x.begin(), x.end());

  double h = 1.0 + (n - 1) * p;
  int hf = std::floor(h);
  double frac = h - hf;
  int i = hf - 1;

  if (i >= n - 1) return x[n - 1];

  return x[i] + frac * (x[i + 1] - x[i]);
}

double trimean_cpp(std::vector<double> x) {
  double q1  = quantile_type7(x, 0.25);
  double med = quantile_type7(x, 0.50);
  double q3  = quantile_type7(x, 0.75);

  return (q1 + 2.0 * med + q3) / 4.0;
}

double median_cpp(std::vector<double> x) {
  return quantile_type7(x, 0.50);
}

double thresholded_mean_cpp(const std::vector<double>& x, double trim) {
  if (x.empty()) return NA_REAL;

  int nnz = 0;
  double sum = 0.0;

  for (double v : x) {
    if (v != 0.0) nnz++;
    sum += v;
  }

  double percent = (double) nnz / x.size();

  if (percent < trim) {
    return 0.0;
  } else {
    return sum / x.size();
  }
}

double truncated_mean_cpp(std::vector<double> x, double trim) {
  int n = x.size();

  if (n == 0) return NA_REAL;
  if (n == 1) return x[0];

  std::sort(x.begin(), x.end());

  int lower = std::floor(n * trim);
  int upper = n - lower;

  if (lower >= upper) return NA_REAL;

  double sum = 0.0;
  int count = 0;

  for (int i = lower; i < upper; i++) {
    sum += x[i];
    count++;
  }

  return sum / count;
}

// [[Rcpp::export]]
NumericMatrix fast_group_summary_cpp(
    NumericMatrix data,
    IntegerVector group,
    int n_groups,
    std::string method = "triMean",
    double trim = 0.1
) {
  int n_genes = data.nrow();
  int n_cells = data.ncol();

  NumericMatrix out(n_genes, n_groups);

  for (int g = 1; g <= n_groups; g++) {
    std::vector<int> cells;
    cells.reserve(n_cells);

    for (int j = 0; j < n_cells; j++) {
      if (group[j] == g) {
        cells.push_back(j);
      }
    }

    for (int i = 0; i < n_genes; i++) {
      std::vector<double> values;
      values.reserve(cells.size());

      for (int idx = 0; idx < (int)cells.size(); idx++) {
        double val = data(i, cells[idx]);

        if (!NumericVector::is_na(val)) {
          values.push_back(val);
        }
      }

      if (method == "triMean") {
        out(i, g - 1) = trimean_cpp(values);
      } else if (method == "truncatedMean") {
        out(i, g - 1) = truncated_mean_cpp(values, trim);
      } else if (method == "thresholdedMean") {
        out(i, g - 1) = thresholded_mean_cpp(values, trim);
      } else if (method == "median") {
        out(i, g - 1) = median_cpp(values);
      } else {
        stop("Unknown method. Use triMean, truncatedMean, thresholdedMean, median, or mean.");
      }
    }
  }

  return out;
}