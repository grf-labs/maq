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

#include <queue>
#include <vector>

#include "compute_path.h"
#include "convex_hull.h"

namespace maq {

struct QueueElement {
  QueueElement(size_t sample, size_t arm, double priority) :
    sample(sample), arm(arm), priority(priority) {}

  size_t sample;
  size_t arm;
  double priority;
};

// break ties by sample
// TODO actually break tie by cost THEN sample?!
// EXACT ties??
bool operator <(const QueueElement& lhs, const QueueElement& rhs) {
  return lhs.priority < rhs.priority
    || (lhs.priority == rhs.priority && lhs.sample > rhs.sample);
}

bool equal_doubles(double first, double second, double epsilon) {
  return std::abs(first - second) < epsilon;
}

solution_path compute_path(const std::vector<size_t>& samples,
                           const std::vector<std::vector<size_t>>& R,
                           const Data& data,
                           double budget,
                           bool bootstrap) {
  std::vector<std::vector<double>> spend_gain(2);
  std::vector<std::vector<size_t>> i_k_path(2);

  std::vector<size_t> active_set(data.num_rows, 0); // active R entry offset by one.
  // std::unordered_map<size_t, size_t> active_set; // slower

  // Initialize PQ with initial enrollment
  std::priority_queue<QueueElement> pqueue;
  for (auto sample : samples) {
    if (!R[sample].empty()) {
      size_t arm = R[sample][0];
      double priority = data.get_reward(sample, arm) / data.get_cost(sample, arm);
      pqueue.emplace(sample, arm, priority);
    }
  }

  double spend = 0;
  double gain = 0;
  while (pqueue.size() > 0 && spend < budget) {
    auto top = pqueue.top();
    pqueue.pop();
    // TODO exact ties

    // assigned before?
    if (active_set[top.sample] > 0) {
      size_t active = active_set[top.sample] - 1;
      size_t active_arm = R[top.sample][active];
      spend -= data.get_cost(top.sample, active_arm);
      gain -= data.get_reward(top.sample, active_arm);
    }

    // assign
    double weight = data.get_weight(top.sample);
    if (bootstrap) {
      // "0-2" bootstrap: half-sample with weight 2.
      weight *= 2;
    }
    double cost = weight * data.get_cost(top.sample, top.arm);
    double reward = weight * data.get_reward(top.sample, top.arm);
    if (spend + cost <= budget || equal_doubles(spend + cost, budget, 1e-16)) {
      spend += cost;
      gain += reward;
      spend_gain[0].push_back(spend);
      spend_gain[1].push_back(gain);
      if (!bootstrap) {
        i_k_path[0].push_back(top.sample);
        i_k_path[1].push_back(top.arm);
      }
      active_set[top.sample]++;
    } else {
      // TODO split decision
      break;
    }

    // upgrade available?
    size_t next_entry = active_set[top.sample];
    if (R[top.sample].size() > next_entry) {
      size_t upgrade = R[top.sample][next_entry];
      double cost_upgrade = data.get_cost(top.sample, upgrade);
      double reward_upgrade = data.get_reward(top.sample, upgrade);
      double priority = (reward_upgrade - reward) / (cost_upgrade - cost);
      pqueue.emplace(top.sample, upgrade, priority);
    }
  }

  // return std::make_pair(std::move(spend_gain), std::move(i_k_path));
  return std::make_pair(spend_gain, i_k_path);
}

} // namespace maq
