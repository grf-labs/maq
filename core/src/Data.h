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

// Data wrapper for column or row major storage
// Rewards can be any real number
// Costs should be > 0
// Weights should be > 0
class Data {
public:
  Data(const double* data_reward,
       const double* data_cost,
       const double* data_weight,
       const int* data_tie_breaker,
       size_t num_rows,
       size_t num_cols,
       bool col_major) :
      num_rows(num_rows),
      num_cols(num_cols),
      data_reward(data_reward),
      data_cost(data_cost),
      data_weight(data_weight),
      data_tie_breaker(data_tie_breaker),
      col_major(col_major) {
    this->has_weight = data_weight == nullptr ? false : true;
    this->has_tie_breaker = data_tie_breaker == nullptr ? false : true;
  }

  double get_reward(size_t row, size_t col) const {
    return data_reward[index(row, col)] * get_weight(row);
  }

  double get_cost(size_t row, size_t col) const {
    return data_cost[index(row, col)] * get_weight(row);
  }

  int get_tie_breaker(size_t row) const {
    if (!has_tie_breaker) {
      return row;
    } else {
      return data_tie_breaker[row];
    }
  }

  size_t num_rows;
  size_t num_cols;

private:
  size_t index(size_t row, size_t col) const {
    if (col_major) {
      return col * num_rows + row;
    } else {
      return row * num_cols + col;
    }
  }

  double get_weight(size_t row) const {
    if (!has_weight) {
      return 1.0 / num_rows;
    } else {
      return data_weight[row];
    }
  }

  const double* data_reward;
  const double* data_cost;
  const double* data_weight;
  const int* data_tie_breaker;
  bool has_weight;
  bool has_tie_breaker;
  bool col_major;
};

} // namespace maq

#endif // MAQ_DATA_H
