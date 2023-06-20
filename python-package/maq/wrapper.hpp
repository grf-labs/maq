#include "MAQ.h"

using namespace maq;

std::pair<solution_path, std::vector<std::vector<double>>> fit(
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

  return run<Storage::RowMajor>(
    data_reward,
    data_reward_scores,
    data_cost,
    num_rows,
    num_cols,
    cost_matrix,
    data_weight,
    data_tie_breaker,
    clusters,
    budget,
    target_with_covariates,
    paired_inference,
    num_bootstrap,
    num_threads,
    random_seed
  );
}
