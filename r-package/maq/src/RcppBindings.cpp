#include <Rcpp.h>

#include "MAQ.h"

using namespace maq;

// [[Rcpp::export]]
Rcpp::List solver_rcpp(const Rcpp::NumericMatrix& reward,
                       const Rcpp::NumericMatrix& reward_scores,
                       const Rcpp::NumericMatrix& cost,
                       const Rcpp::NumericVector& sample_weights,
                       const Rcpp::IntegerVector& tie_breaker,
                       const Rcpp::IntegerVector& clusters,
                       double budget,
                       unsigned int num_bootstrap,
                       unsigned int num_threads,
                       unsigned int seed) {
  size_t num_rows = reward.rows();
  size_t num_cols = reward.cols();
  const double* weights_ptr = nullptr;
  if (sample_weights.size() > 0) {
    weights_ptr = sample_weights.begin();
  }
  const int* tie_breaker_ptr = nullptr;
  if (tie_breaker.size() > 0) {
    tie_breaker_ptr = tie_breaker.begin();
  }
  const int* clusters_ptr = nullptr;
  if (clusters.size() > 0) {
    clusters_ptr = clusters.begin();
  }

  Data data(reward.begin(), reward_scores.begin(), cost.begin(), weights_ptr, tie_breaker_ptr, clusters_ptr, num_rows, num_cols, true);
  MAQOptions options(budget, num_bootstrap, num_threads, seed);
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
  Data data(reward.begin(), reward.begin(), cost.begin(), nullptr, nullptr, nullptr, num_rows, num_cols, true);

  return Rcpp::List::create(convex_hull(data));
}
