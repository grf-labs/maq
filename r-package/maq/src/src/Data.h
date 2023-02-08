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

namespace maq {

// Data wrapper for column major storage
class Data {
public:
  Data(const double* data_reward,
       const double* data_cost,
       const int* data_order,
       const double* data_weight,
       size_t num_rows,
       size_t num_cols) :
      num_rows(num_rows),
      num_cols(num_cols),
      data_reward(data_reward),
      data_cost(data_cost),
      data_order(data_order),
      data_weight(data_weight) {
    if (data_weight == nullptr) {
      this->has_weight = false;
    } else {
      this->has_weight = true;
    }
    // double weight_sum = 0;
    // for (size_t i = 0; i < num_rows; i++) {
    //   weight_sum += get_weight(i);
    // }
    // this->weight_sum = weight_sum;
  }

  // Rewards can be any real number
  double get_reward(size_t row, size_t col) const {
    return data_reward[col * num_rows + row];
  }

  // Costs should be > 0
  double get_cost(size_t row, size_t col) const {
    return data_cost[col * num_rows + row];
  }

  // Return the sample index of points in order of increasing cost.
  // data_order should be a num_rows * num_cols array of the row-wise
  // sort order of data_cost (0-indexed).
  size_t get_order(size_t row, size_t col) const {
    return data_order[col * num_rows + row];
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
  // double weight_sum;

private:
  const double* data_reward;
  const double* data_cost;
  const int* data_order;
  const double* data_weight;
  bool has_weight;
};

} // namespace maq

#endif // MAQ_DATA_H
