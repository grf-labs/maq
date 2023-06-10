#ifndef MAQ_COMPLETEDATA_H
#define MAQ_COMPLETEDATA_H

#include <cstddef>
#include <vector>

#include "Data.hpp"

namespace maq {

template <class T>
class CompleteData{
public:
  CompleteData(const T& data) :
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
  const T& data;
};

} // namespace maq

#endif // MAQ_COMPLETEDATA_H
