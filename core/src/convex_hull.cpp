/*-------------------------------------------------------------------------------
  This file is part of Multi-Action Qini (maq).

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

#include <algorithm>
#include <numeric>
#include <vector>

#include "convex_hull.h"

/*
Find the upper left convex hull on the (cost, reward) plane for each sample. This takes
O(num_rows * num_arms * log(num_arms)) time using the Graham scan for finding a convex hull,
with the the "angle calculation" replaced by a variant of the "LP dominance" criteria
(11.6) and (11.7) in Kellerer et al. (2004, Chapter 11).

Consider 3 points j, k, l ordered according to cost j < k < l.
This algorithm maintains a stack with the first two points on top, then iteratively
proceeds by checking if point k should be replaced by or augmented by point l.
*/

namespace maq {

inline bool is_dominated(const std::vector<size_t>& Ri,
                         size_t point_l,
                         size_t sample,
                         const Data& data) {
  if (Ri.size() < 1) {
    return false;
  }
  // a dummy (0, 0) origin point emulating a zero-cost control.
  double cost_j = 0;
  double reward_j = 0;
  if (Ri.size() >= 2) {
    size_t point_j = Ri[Ri.size() - 2]; // next to top
    cost_j = data.get_cost(sample, point_j);
    reward_j = data.get_reward(sample, point_j);
  }
  size_t point_k = Ri[Ri.size() - 1]; // top
  double cost_k = data.get_cost(sample, point_k);
  double reward_k = data.get_reward(sample, point_k);
  if (reward_k <= 0) {
    return true;
  }
  double cost_l = data.get_cost(sample, point_l);
  double reward_l = data.get_reward(sample, point_l);

  // C++: a/0 = Inf if a > 0, a/0 = -Inf if a <0, and 0/0 = NaN (all logical operators on NaN evaluate to false)
  return (reward_l - reward_k) / (cost_l - cost_k) > (reward_k - reward_j) / (cost_k - cost_j);
}

std::vector<std::vector<size_t>> convex_hull(const Data& data) {
  std::vector<std::vector<size_t>> R(data.num_rows);
  std::vector<size_t> ordered_arms(data.num_cols);
  std::iota(ordered_arms.begin(), ordered_arms.end(), 0); // fill with 0, ..., K - 1

  for (size_t sample = 0; sample < data.num_rows; sample++) {
    std::vector<size_t>& Ri = R[sample];

    // Get sort order by increasing cost
    std::stable_sort(ordered_arms.begin(), ordered_arms.end(), [&](const size_t lhs, const size_t rhs) {
      return data.get_cost(sample, lhs) < data.get_cost(sample, rhs);
    });
    // Push first positive reward point onto stack
    size_t start = 0;
    while (start < data.num_cols && data.get_reward(sample, ordered_arms[start]) <= 0) {
      start++;
    }
    if (start == data.num_cols) {
      continue;
    }
    size_t first = ordered_arms[start];
    Ri.push_back(first);
    for (size_t l = start + 1; l < data.num_cols; l++) {
      size_t point_l = ordered_arms[l];
      while (is_dominated(Ri, point_l, sample, data)) {
        Ri.pop_back(); // remove point_k
      }
      double reward_l = data.get_reward(sample, point_l);
      if (reward_l > 0) {
        if (Ri.size() < 1 || reward_l > data.get_reward(sample, Ri.back())) {
          Ri.push_back(point_l);
        }
      }
    }
  }

  return R;
}

} // namespace maq
