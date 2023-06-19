#ifndef MAQ_H
#define MAQ_H

// Multi-Armed Qini (maq).
// https://github.com/grf-labs/maq
// Distributed under the MIT License.

#include "Data.hpp"
#include "SolverOptions.hpp"
#include "Solver.hpp"

namespace maq {

template <class T>
Solver<T> make_solver(const T& data, const SolverOptions& options) {return Solver<T>(data, options);}

template <Storage storage>
std::pair<solution_path, std::vector<std::vector<double>>> run(
    const double* data_reward,
    const double* data_reward_scores,
    const double* data_cost,
    size_t num_rows,
    size_t num_cols,
    bool cost_matrix,
    const double* data_weight,
    const int* data_tie_breaker,
    const int* clusters,
    double budget,
    bool target_with_covariates,
    bool paired_inference,
    unsigned int num_bootstrap,
    unsigned int num_threads,
    unsigned int random_seed) {
  SolverOptions options(budget, target_with_covariates, paired_inference, num_bootstrap, num_threads, random_seed);

  if (data_weight == nullptr && data_tie_breaker == nullptr) {
    Data<storage, SampleWeights::Default, TieBreaker::Default> data(
      data_reward, data_reward_scores, data_cost, num_rows, num_cols,
      data_weight, data_tie_breaker, clusters);
    auto maq = make_solver(data, options);
    return maq.fit();
  } else if (data_weight != nullptr && data_tie_breaker == nullptr) {
    Data<storage, SampleWeights::Provided, TieBreaker::Default> data(
      data_reward, data_reward_scores, data_cost, num_rows, num_cols,
      data_weight, data_tie_breaker, clusters);
    auto maq = make_solver(data, options);
    return maq.fit();
  } else if (data_weight != nullptr && data_tie_breaker != nullptr) {
    Data<storage, SampleWeights::Provided, TieBreaker::Provided> data(
      data_reward, data_reward_scores, data_cost, num_rows, num_cols,
      data_weight, data_tie_breaker, clusters);
    auto maq = make_solver(data, options);
    return maq.fit();
  } else {
    Data<storage, SampleWeights::Default, TieBreaker::Provided> data(
      data_reward, data_reward_scores, data_cost, num_rows, num_cols,
      data_weight, data_tie_breaker, clusters);
    auto maq = make_solver(data, options);
    return maq.fit();
  }
}


} // namespace maq

#endif // MAQ_H
