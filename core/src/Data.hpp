#ifndef MAQ_DATA_H
#define MAQ_DATA_H

#include <cstddef>
#include <vector>

namespace maq {

enum class Storage {ColMajor, RowMajor};
enum class SampleWeights {Default, Provided};
enum class TieBreaker {Default, Provided};

/**
 * Read-only data wrapper for column or row major storage.
 *
 * Costs should be > 0
 * Weights should be > 0 and sum to 1
 * Clusters, if present, should be labeled as consecutive integers 0, ..., num_clusters
 *
 */
template <Storage storage, SampleWeights sample_weights, TieBreaker tie_breaker>
class Data {
public:
  Data(const double* data_reward,
       const double* data_reward_scores,
       const double* data_cost,
       size_t num_rows,
       size_t num_cols,
       const double* data_weight = nullptr,
       const int* data_tie_breaker = nullptr,
       const int* clusters = nullptr) :
    num_rows(num_rows),
    num_cols(num_cols),
    data_reward(data_reward),
    data_reward_scores(data_reward_scores),
    data_cost(data_cost),
    data_weight(data_weight),
    data_tie_breaker(data_tie_breaker) {

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

  double get_reward(size_t row, size_t col) const {
     return data_reward[index(row, col)] * get_weight(row);
  }

  double get_reward_scores(size_t row, size_t col) const {
    return data_reward_scores[index(row, col)] * get_weight(row);
  }

  double get_cost(size_t row, size_t col) const {
    return data_cost[index(row, col)] * get_weight(row);
  }

  int get_tie_breaker(size_t row) const {
    if (tie_breaker == TieBreaker::Default) {
    // if (data_tie_breaker == nullptr) {
      return row;
    } else {
      return data_tie_breaker[row];
    }
  }

  size_t num_rows;
  size_t num_cols;
  std::vector<std::vector<size_t>> samples_by_cluster;

private:
  double get_weight(size_t row) const {
    if (sample_weights == SampleWeights::Default) {
    // if (data_weight == nullptr) {
      return 1.0 / num_rows;
      } else {
      return data_weight[row];
    }
  }

  size_t index(size_t row, size_t col) const {
    if (storage == Storage::ColMajor) {
    // if (row >= 0) {
      return col * num_rows + row;
    } else {
      return row * num_cols + col;
    }
  }


  const double* data_reward;
  const double* data_reward_scores;
  const double* data_cost;
  const double* data_weight;
  const int* data_tie_breaker;
};

} // namespace maq

#endif // MAQ_DATA_H
