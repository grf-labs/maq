#include <cmath>
#include <queue>
#include <vector>

#include "compute_path.h"
#include "convex_hull.h"

namespace maq {

struct QueueElement {
  QueueElement(size_t sample, size_t arm, int tie_breaker, double priority) :
    sample(sample), arm(arm), tie_breaker(tie_breaker), priority(priority) {}

  size_t sample;
  size_t arm;
  int tie_breaker;
  double priority;
};

bool operator <(const QueueElement& lhs, const QueueElement& rhs) {
  return lhs.priority < rhs.priority
    || (lhs.priority == rhs.priority && lhs.tie_breaker > rhs.tie_breaker);
}

solution_path compute_path(const std::vector<size_t>& samples,
                           const std::vector<std::vector<size_t>>& R,
                           const Data& data,
                           double budget,
                           bool bootstrap) {
  std::vector<std::vector<double>> spend_gain(3); // 3rd entry: SEs
  std::vector<std::vector<size_t>> i_k_path(3); // 3rd entry: complete path?
  std::vector<size_t> active_set(data.num_rows, 0); // active R entry offset by one (vec faster than hash table)

  // Initialize PQ with initial enrollment
  std::priority_queue<QueueElement> pqueue;
  for (auto sample : samples) {
    if (!R[sample].empty()) {
      size_t arm = R[sample][0];
      int tie_breaker = data.get_tie_breaker(sample);
      double priority = data.get_reward(sample, arm) / data.get_cost(sample, arm);
      pqueue.emplace(sample, arm, tie_breaker, priority);
    }
  }

  double spend = 0;
  double gain = 0;
  double bs_weight = bootstrap ? 2 : 1; // "0-2" bootstrap: half-sample with weight 2.
  while (pqueue.size() > 0 && spend < budget) {
    auto top = pqueue.top();
    pqueue.pop();

    // assigned before?
    if (active_set[top.sample] > 0) {
      size_t active = active_set[top.sample] - 1;
      size_t active_arm = R[top.sample][active];
      spend -= bs_weight * data.get_cost(top.sample, active_arm);
      gain -= bs_weight * data.get_reward_scores(top.sample, active_arm);
    }

    // assign
    double cost = data.get_cost(top.sample, top.arm);
    double reward = data.get_reward(top.sample, top.arm);

    spend += bs_weight * cost;
    gain += bs_weight * data.get_reward_scores(top.sample, top.arm);
    spend_gain[0].push_back(spend);
    spend_gain[1].push_back(gain);
    if (!bootstrap) {
      i_k_path[0].push_back(top.sample);
      i_k_path[1].push_back(top.arm);
    }
    active_set[top.sample]++;

    // upgrade available?
    size_t next_entry = active_set[top.sample];
    if (R[top.sample].size() > next_entry) {
      size_t upgrade = R[top.sample][next_entry];
      double cost_upgrade = data.get_cost(top.sample, upgrade);
      double reward_upgrade = data.get_reward(top.sample, upgrade);
      double priority = (reward_upgrade - reward) / (cost_upgrade - cost);
      pqueue.emplace(top.sample, upgrade, top.tie_breaker, priority);
    }

    // have we reached maximum spend? if so stop at nearest integer solution (rounded up)
    if (spend >= budget) {
      break;
    }
  }

  if (!bootstrap) {
    i_k_path[2].push_back(pqueue.size() > 0 ? 1 : 0);
  }

  return std::make_pair(std::move(spend_gain), std::move(i_k_path));
}

} // namespace maq
