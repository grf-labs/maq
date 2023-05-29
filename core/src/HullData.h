#ifndef MAQ_HULLDATA_H
#define MAQ_HULLDATA_H

#include <cstddef>
#include <vector>

#include "Data.h"

namespace maq {

class HullData {
public:
  virtual double get_reward(size_t row, size_t col) const = 0;

  virtual double get_cost(size_t row, size_t col) const = 0;

  virtual size_t get_num_rows() const = 0;

  virtual size_t get_num_cols() const = 0;

  virtual double get_reward_scores(size_t row, size_t col) const = 0;
};

class CompleteData : public HullData {
public:
  CompleteData(const Data& data) :
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
  const Data& data;
};

class MeanData : public HullData {
public:
  MeanData(const Data& data, const std::vector<size_t>& samples) {
    std::vector<double> reward(data.num_cols);
    std::vector<double> reward_scores(data.num_cols);
    std::vector<double> cost(data.num_cols);
    for (auto sample : samples) {
      for (size_t col = 0; col < data.num_cols; col++) {
        reward[col] += data.get_reward(sample, col);
        reward_scores[col] += data.get_reward_scores(sample, col);
        cost[col] += data.get_cost(sample, col);
      }
    }
    this->num_rows = data.num_rows;
    this->reward = reward;
    this->reward_scores = reward_scores;
    this->cost = cost;
  }

  double get_reward(size_t row, size_t col) const {
    return reward[col] / num_rows;
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

  double get_reward_scores(size_t row, size_t col) const {
    return reward_scores[col] / num_rows;
  }

private:
  size_t num_rows;
  std::vector<double> reward;
  std::vector<double> reward_scores;
  std::vector<double> cost;
};

} // namespace maq

#endif // MAQ_HULLDATA_H
