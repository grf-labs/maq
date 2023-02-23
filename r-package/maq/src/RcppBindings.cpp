/*-------------------------------------------------------------------------------
  This file is part of Multi-Action QINI (maq).

  maq is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  maq is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with maq. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#include <Rcpp.h>

#include "MAQ.h"

using namespace maq;

// [[Rcpp::export]]
Rcpp::List solver_rcpp(const Rcpp::NumericMatrix& reward,
                       const Rcpp::NumericMatrix& cost,
                       const Rcpp::NumericVector& sample_weights,
                       const std::vector<size_t>& clusters,
                       uint samples_per_cluster,
                       double budget,
                       size_t num_bootstrap,
                       uint num_threads,
                       uint seed) {
  size_t num_rows = reward.rows();
  size_t num_cols = reward.cols();
  const double* weights_ptr = nullptr;
  if (sample_weights.size() > 0) {
    weights_ptr = sample_weights.begin();
  }

  Data data(reward.begin(), cost.begin(), weights_ptr, num_rows, num_cols);
  MAQOptions options(budget, num_bootstrap, clusters, samples_per_cluster, num_threads, seed);
  MAQ maq(data, options);

  auto ret = maq.fit();

  Rcpp::List res;
  res.push_back(ret.first[0], "spend");
  res.push_back(ret.first[1], "gain");
  res.push_back(ret.first[2], "std.err");
  res.push_back(ret.second[0], "ipath");
  res.push_back(ret.second[1], "kpath");
  res.push_back(ret.second[2][0] > 0 ? false : true, "complete.path");

  return res;
}

// this function is only wrapped for testing purposes.
#include "convex_hull.h"
// [[Rcpp::export]]
Rcpp::List convex_hull_rcpp(const Rcpp::NumericMatrix& reward,
                            const Rcpp::NumericMatrix& cost) {
  size_t num_rows = reward.rows();
  size_t num_cols = reward.cols();
  Data data(reward.begin(), cost.begin(), nullptr, num_rows, num_cols);

  return Rcpp::List::create(convex_hull(data));
}
