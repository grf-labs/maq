/*-------------------------------------------------------------------------------
  This file is part of Multi-Action QINI (maq).

  maq is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  maq is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with maq. If not, see <http://www.gnu.org/licenses/>.
 #-------------------------------------------------------------------------------*/

#ifndef MAQ_DATA_H
#define MAQ_DATA_H

#include <cstddef>

namespace maq {

// Data wrapper for column major storage
// Rewards can be any real number
// Costs should be > 0
class Data {
public:
  Data(const double* data_reward,
       const double* data_cost,
       const double* data_weight,
       size_t num_rows,
       size_t num_cols) :
      num_rows(num_rows),
      num_cols(num_cols),
      data_reward(data_reward),
      data_cost(data_cost),
      data_weight(data_weight) {
    if (data_weight == nullptr) {
      this->has_weight = false;
    } else {
      this->has_weight = true;
    }
    double weight_sum = 0;
    for (size_t i = 0; i < num_rows; i++) {
      weight_sum += get_weight(i);
    }
    this->weight_sum = weight_sum;
  }

  double get_reward(size_t row, size_t col) const {
    return data_reward[col * num_rows + row] / weight_sum;
  }

  double get_cost(size_t row, size_t col) const {
    return data_cost[col * num_rows + row] / weight_sum;;
  }

  double get_weight(size_t row) const {
    if (!has_weight) {
      return 1.0;
    } else {
      return data_weight[row];
    }
  }

  size_t num_rows;
  size_t num_cols;

private:
  const double* data_reward;
  const double* data_cost;
  const double* data_weight;
  bool has_weight;
  double weight_sum;
};

} // namespace maq

#endif // MAQ_DATA_H
