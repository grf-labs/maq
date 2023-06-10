#ifndef MAQ_COMPLETEDATA_H
#define MAQ_COMPLETEDATA_H

#include <cstddef>
#include <vector>

#include "Data.h"

namespace maq {

template <Storage storage, SampleWeights sample_weights, TieBreaker tie_breaker>
class CompleteData{
public:
  CompleteData(const Data<storage, sample_weights, tie_breaker>& data) :
    data(data) {}

  double get_reward(size_t row, size_t col) const {
    return data.get_reward(row, col);
  }

  double get_cost(size_t row, size_t col) const {
    return data.get_cost(row, col);
  }

  size_t get_num_rows() const {
    return data.num_rows;
  }

  size_t get_num_cols() const {
    return data.num_cols;
  }

  double get_reward_scores(size_t row, size_t col) const {
    return data.get_reward_scores(row, col);
  }

private:
  const Data<storage, sample_weights, tie_breaker>& data;
};

} // namespace maq

#endif // MAQ_COMPLETEDATA_H
