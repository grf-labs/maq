#include "MAQ.h"

using namespace maq;

solution_path run(const double* data_reward,
                  const double* data_reward_scores,
                  const double* data_cost,
                  size_t num_rows,
                  size_t num_cols,
                  double budget,
                  bool target_with_covariates,
                  unsigned int num_bootstrap,
                  unsigned int num_threads,
                  unsigned int random_seed) {
  SolverOptions options(budget, target_with_covariates, false, num_bootstrap, num_threads, random_seed);
  Data<Storage::RowMajor> data(
    data_reward, data_reward_scores, data_cost, num_rows, num_cols);
  auto maq = make_solver(data, options);

  return maq.fit().first;
}
