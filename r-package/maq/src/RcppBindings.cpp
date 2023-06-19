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
                       bool target_with_covariates,
                       bool paired_inference,
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
  auto ret = run<Storage::ColMajor>(
    reward.begin(),
    reward_scores.begin(),
    cost.begin(),
    num_rows,
    num_cols,
    true,
    weights_ptr,
    tie_breaker_ptr,
    clusters_ptr,
    budget,
    target_with_covariates,
    paired_inference,
    num_bootstrap,
    num_threads,
    seed
  );
  auto path = ret.first;

  Rcpp::List res;
  res.push_back(path.first[0], "spend");
  res.push_back(path.first[1], "gain");
  res.push_back(path.first[2], "std.err");
  res.push_back(path.second[0], "ipath");
  res.push_back(path.second[1], "kpath");
  res.push_back(path.second[2][0] > 0 ? true : false, "complete.path");
  res.push_back(ret.second, "gain.bs");

  return res;
}

// this function is only wrapped for testing purposes.
// [[Rcpp::export]]
Rcpp::List convex_hull_rcpp(const Rcpp::NumericMatrix& reward,
                            const Rcpp::NumericMatrix& cost) {
  size_t num_rows = reward.rows();
  size_t num_cols = reward.cols();
  Data<Storage::ColMajor, SampleWeights::Default, TieBreaker::Default> data(
    reward.begin(), reward.begin(), cost.begin(), num_rows, num_cols);

  return Rcpp::List::create(convex_hull(data));
}
