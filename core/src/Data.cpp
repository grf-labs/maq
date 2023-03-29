#include "Data.h"

namespace maq {

Data::Data(const double* data_reward,
       const double* data_reward_scores,
       const double* data_cost,
       const double* data_weight,
       const int* data_tie_breaker,
       const int* clusters,
       size_t num_rows,
       size_t num_cols,
       bool col_major) :
      num_rows(num_rows),
      num_cols(num_cols),
      data_reward(data_reward),
      data_reward_scores(data_reward_scores),
      data_cost(data_cost),
      data_weight(data_weight),
      data_tie_breaker(data_tie_breaker),
      col_major(col_major) {

  this->has_weight = data_weight == nullptr ? false : true;
  this->has_tie_breaker = data_tie_breaker == nullptr ? false : true;

  // If clusters are present, then fill samples_by_cluster with samples belonging to each cluster.
  if (clusters != nullptr) {
    for (size_t sample = 0; sample < num_rows; sample++) {
      size_t cluster_id = clusters[sample];
      if (cluster_id + 1 > samples_by_cluster.size()) {
        samples_by_cluster.resize(cluster_id + 1);
      }
      samples_by_cluster[cluster_id].push_back(sample);
    }
  }
}

} // namespace maq
