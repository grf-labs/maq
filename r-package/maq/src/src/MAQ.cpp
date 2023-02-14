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

#include <future>
#include <random>
#include <vector>

#include "MAQ.h"
#include "compute_path.h"
#include "convex_hull.h"
#include "sampling/SamplingOptions.h"
#include "sampling/RandomSampler.h"

// #include<iostream>

namespace maq {

MAQ::MAQ(const Data& data,
         const MAQOptions& options) :
    data(data),
    options(options) {
  this->sampling_options = grf::SamplingOptions(options.samples_per_cluster, options.clusters);
}

solution_path MAQ::fit() {
  std::vector<size_t> samples;
  samples.reserve(data.num_rows);
  for (size_t sample = 0; sample < data.num_rows; sample++) {
    samples.push_back(sample);
  }

  auto R = convex_hull(data);
  auto path_hat = compute_path(samples, R, data, options.budget, false);
  auto gain_interp = fit_paths(path_hat, R);
  compute_std_err(path_hat, gain_interp);

  return path_hat;
}

std::vector<std::vector<double>> MAQ::fit_paths(const solution_path& path_hat,
                                                const std::vector<std::vector<size_t>>& R) {
  std::vector<uint> thread_ranges;
  split_sequence(thread_ranges, 0, static_cast<uint>(options.num_bootstrap - 1), options.num_threads);

  std::vector<std::future<std::vector<std::vector<double>>>> futures;
  futures.reserve(thread_ranges.size());

  std::vector<std::vector<double>> predictions;
  predictions.reserve(options.num_bootstrap);

  for (uint i = 0; i < thread_ranges.size() - 1; ++i) {
    size_t start_index = thread_ranges[i];
    size_t num_replicates_batch = thread_ranges[i + 1] - start_index;

    futures.push_back(std::async(std::launch::async,
                                 &MAQ::fit_paths_batch,
                                 this,
                                 start_index,
                                 num_replicates_batch,
                                 std::ref(path_hat),
                                 std::ref(R)));
  }

  for (auto& future : futures) {
    auto thread_predictions = future.get();
    predictions.insert(predictions.end(),
                       std::make_move_iterator(thread_predictions.begin()),
                       std::make_move_iterator(thread_predictions.end()));
  }

  return predictions;
}

std::vector<std::vector<double>> MAQ::fit_paths_batch(size_t start,
                                                      size_t num_replicates,
                                                      const solution_path& path_hat,
                                                      const std::vector<std::vector<size_t>>& R) {
  std::mt19937_64 random_number_generator(options.random_seed + start);
  nonstd::uniform_int_distribution<uint> udist;

  std::vector<std::vector<double>> predictions;
  predictions.reserve(num_replicates);

  for (size_t b = 0; b < num_replicates; b++) {
    uint bs_seed = udist(random_number_generator);
    grf::RandomSampler sampler(bs_seed, sampling_options);
    std::vector<size_t> clusters;
    std::vector<size_t> samples;
    sampler.sample_clusters(data.num_rows, 0.5, clusters);
    sampler.sample_from_clusters(clusters, samples);
    // sampler.sample_clusters(data.num_rows, 0.5, samples);
    // sampler.sample_from_clusters(samples, samples);

    auto path_b = compute_path(samples, R, data, options.budget, true);
    auto gain_b = interpolate_path(path_hat, path_b);
    predictions.push_back(std::move(gain_b));
  }

  return predictions;
}

void MAQ::compute_std_err(solution_path& path_hat, const std::vector<std::vector<double>>& gain_interp) {
  size_t grid_len = path_hat.first[0].size();
  std::vector<double>& std_err = path_hat.first[2];
  std_err.resize(grid_len);
  if (gain_interp.size() < 2) {
    return;
  }

  for (size_t i = 0; i < grid_len; i++) {
    // Use Welford's algorithm to get numerically stable variance estimates in one pass.
    double Mprev;
    double M;
    double Sprev = -1;
    double S;
    size_t n = 0;
    for (size_t b = 0; b < gain_interp.size(); b++) {
      if (gain_interp[b].size() < 1) {
        continue;
      }
      double val = gain_interp[b][i];
      if (val == -1) {
        continue;
      }
      n++;
      if (Sprev == -1) {
        Mprev = val;
        Sprev = 0;
        continue;
      }
      M = Mprev + (val - Mprev) / n;
      S = Sprev + (val - Mprev) * (val - M);

      Mprev = M;
      Sprev = S;
    }
    if (n >= 2) {
      std_err[i] = sqrt(S / (n - 1));
    } else {
      std_err[i] = 0; // define these to be zero.
    }
  }

}

std::vector<double> MAQ::interpolate_path(const solution_path& path_hat, const solution_path& path_hat_b) {
  // interpolate bootstrapped gain on \hat path's (monotonically increasing) spend grid.
  const std::vector<double>& grid = path_hat.first[0];
  std::vector<double> interp;

  const std::vector<double>& grid_b = path_hat_b.first[0];
  const std::vector<double>& gain_b = path_hat_b.first[1];
  if (grid_b.size() < 1) {
    return interp;
  }
  interp.resize(grid.size());

  // initialize the interpolation interval
  size_t left = 0;
  size_t right = grid_b.size() < 2 ? 0 : 1;
  for (size_t i = 0; i < grid.size(); i++) {
    double val = grid[i];
    // out of left range?
    if (val < grid_b[left]) {
      interp[i] = -1;
      continue;
    }
    // update active interval?
    while (right + 2 <= grid_b.size() && grid_b[left + 1] <= val) {
      left++;
      right++;
    }
    // out of right range?
    if (val >= grid_b[right]) {
      interp[i] = gain_b[right];
      continue;
    }
    interp[i] = gain_b[left] + (gain_b[right] - gain_b[left]) *
                  (val - grid_b[left]) / (grid_b[right] - grid_b[left]);
  }

  return interp;
}

void MAQ::split_sequence(std::vector<uint>& result, uint start, uint end, uint num_parts) {
  result.reserve(num_parts + 1);

  // Return range if only 1 part
  if (num_parts == 1) {
    result.push_back(start);
    result.push_back(end + 1);
    return;
  }

  // Return vector from start to end+1 if more parts than elements
  if (num_parts > end - start + 1) {
    for (uint i = start; i <= end + 1; ++i) {
      result.push_back(i);
    }
    return;
  }

  uint length = (end - start + 1);
  uint part_length_short = length / num_parts;
  uint part_length_long = (uint) std::ceil(length / ((double) num_parts));
  uint cut_pos = length % num_parts;

  // Add long ranges
  for (uint i = start; i < start + cut_pos * part_length_long; i = i + part_length_long) {
    result.push_back(i);
  }

  // Add short ranges
  for (uint i = start + cut_pos * part_length_long; i <= end + 1; i = i + part_length_short) {
    result.push_back(i);
  }
}

} // namespace maq
