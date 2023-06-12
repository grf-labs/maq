#ifndef MAQ_DATAMEAN_HPP
#define MAQ_DATAMEAN_HPP

#include <cstddef>
#include <vector>

namespace maq {

template <class DataType>
class DataMean {
public:
  DataMean(const DataType& data, const std::vector<size_t>& samples) {
    reward.resize(data.get_num_cols());
    reward_scores.resize(data.get_num_cols());
    cost.resize(data.get_num_cols());
    for (auto sample : samples) {
      for (size_t col = 0; col < data.get_num_cols(); col++) {
        reward[col] += data.get_reward(sample, col);
        reward_scores[col] += data.get_reward_scores(sample, col);
        cost[col] += data.get_cost(sample, col);
      }
    }
    num_rows = data.get_num_rows();
  }

  double get_reward(size_t row, size_t col) const {
    return reward[col] / num_rows;
  }

  double get_reward_scores(size_t row, size_t col) const {
    return reward_scores[col] / num_rows;
  }

  double get_cost(size_t row, size_t col) const {
    return cost[col] / num_rows;
  }

  size_t get_num_rows() const {
    return 1;
  }

  size_t get_num_cols() const {
    return reward.size();
  }

private:
  size_t num_rows;
  std::vector<double> reward;
  std::vector<double> reward_scores;
  std::vector<double> cost;
};

} // namespace maq

#endif // MAQ_DATAMEAN_HPP
