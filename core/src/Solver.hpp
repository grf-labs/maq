#ifndef MAQ_SOLVER_HPP
#define MAQ_SOLVER_HPP

// Multi-Armed Qini (maq).
// https://github.com/grf-labs/maq
// Distributed under the MIT License.

#include <cmath>
#include <future>
#include <random>
#include <vector>

#include "convex_hull.hpp"
#include "compute_path.hpp"
#include "DataMean.hpp"
#include "Sampler.hpp"
#include "SolverOptions.hpp"

namespace maq {

/**
 * Fit a Multi-Armed Qini curve.
 *
 * The solution_path is a pair where the first entry is a set of vectors containing the path
 * {spend, gain, std.err} and the second pair the path of the corresponding allocations {unit
 * index, arm index}. The second entry in fit()’s return pair is a set of vectors with the
 * bootstrapped gain path. Since the spend path in general can be an irregular grid (each
 * unit/arm may have different cost), these paths are interpolated and indexed to the original
 * full-sample spend path. This entry is non-empty if options.paired_inference = true, and
 * enables paired-wise tests with other cost curves that are fit with the same bootstrap seed
 * (thus yielding valid standard errors that account for the correlation of two curves estimated
 * on the same bootstrap sample).
 *
 * The vectors {spend, gain} trace out the path of integer solutions to the underlying LP, and
 * essentially act as a sufficient statistic that gives the solution to any point on the cost
 * curves. For example, the LP solution for a budget  equal to some arbitrary point between
 * spend[j] and spend[j+1] is given by the linear interpolation between gain[j] and gain[j+1].
 * Likewise, the vectors {unit index, arm index} implicitly act as similar compressed statistics
 * that can yield the solution \pi matrix at any budget level. These steps are performed in the
 * wrapper code.
 *
 */
template <class DataType>
class Solver {
  public:
    Solver(const DataType& data,
           const SolverOptions& options) : data(data), options(options) {}

    std::pair<solution_path, std::vector<std::vector<double>>> fit() {
      std::vector<size_t> samples;
      samples.reserve(data.get_num_rows());
      for (size_t sample = 0; sample < data.get_num_rows(); sample++) {
        samples.push_back(sample);
      }

      std::vector<std::vector<size_t>> R;
      solution_path path_hat;
      if (options.target_with_covariates) {
        R = convex_hull(data);
        path_hat = compute_path(samples, R, data, options.budget, false);
      } else {
        auto mean_data = DataMean<DataType>(data, samples);
        R = convex_hull(mean_data);
        path_hat = compute_path(samples, R[0], mean_data, options.budget, false);
      }

      auto gain_bs = fit_paths(path_hat, R);
      compute_std_err(path_hat, gain_bs);

      return std::make_pair(path_hat,
        options.paired_inference ? std::move(gain_bs) : std::vector<std::vector<double>>());
    }

  private:
    std::vector<std::vector<double>> fit_paths(const solution_path& path_hat,
                                               const std::vector<std::vector<size_t>>& R) {
      std::vector<unsigned int> thread_ranges;
      split_sequence(thread_ranges, 0, static_cast<unsigned int>(options.num_bootstrap - 1), options.num_threads);

      std::vector<std::future<std::vector<std::vector<double>>>> futures;
      futures.reserve(thread_ranges.size());

      std::vector<std::vector<double>> predictions;
      predictions.reserve(options.num_bootstrap);

      for (unsigned int i = 0; i < thread_ranges.size() - 1; ++i) {
        size_t start_index = thread_ranges[i];
        size_t num_replicates_batch = thread_ranges[i + 1] - start_index;

        futures.push_back(std::async(std::launch::async,
                                    &Solver::fit_paths_batch,
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

    std::vector<std::vector<double>> fit_paths_batch(size_t start,
                                                     size_t num_replicates,
                                                     const solution_path& path_hat,
                                                     const std::vector<std::vector<size_t>>& R) {
      std::vector<std::vector<double>> predictions;
      predictions.reserve(num_replicates);

      for (size_t b = 0; b < num_replicates; b++) {
        std::vector<size_t> samples = Sampler<DataType>::sample(data, 0.5, options.random_seed + start + b);
        solution_path path_b;
        if (options.target_with_covariates) {
          path_b = compute_path(samples, R, data, options.budget, true);
        } else {
          auto mean_data = DataMean<DataType>(data, samples);
          auto R_mean = convex_hull(mean_data);
          path_b = compute_path(samples, R_mean[0], mean_data, options.budget, true);
        }
        auto gain_b = interpolate_path(path_hat, path_b);
        predictions.push_back(std::move(gain_b));
      }

      return predictions;
    }

    std::vector<double> interpolate_path(const solution_path& path_hat,
                                         const solution_path& path_hat_b) {
      // interpolate bootstrapped gain on \hat path's (monotonically increasing) spend grid.
      const std::vector<double>& grid = path_hat.first[0];
      std::vector<double> gain_b_interp;

      const std::vector<double>& grid_b = path_hat_b.first[0];
      const std::vector<double>& gain_b = path_hat_b.first[1];
      if (grid_b.size() < 1) {
        return gain_b_interp;
      }
      gain_b_interp.resize(grid.size());

      // initialize the interpolation interval
      size_t left = 0;
      size_t right = grid_b.size() < 2 ? 0 : 1;
      for (size_t i = 0; i < grid.size(); i++) {
        double val = grid[i];
        // out of left range?
        if (val < grid_b[left]) {
          gain_b_interp[i] = NAN;
          continue;
        }
        // update active interval?
        while (right + 2 <= grid_b.size() && grid_b[left + 1] <= val) {
          left++;
          right++;
        }
        // out of right range?
        if (val >= grid_b[right]) {
          gain_b_interp[i] = gain_b[right];
          continue;
        }
        gain_b_interp[i] = gain_b[left] + (gain_b[right] - gain_b[left]) *
                            (val - grid_b[left]) / (grid_b[right] - grid_b[left]);
      }

      return gain_b_interp;
    }

    void compute_std_err(solution_path& path_hat,
                         const std::vector<std::vector<double>>& gain_bs) {
      size_t grid_len = path_hat.first[0].size();
      std::vector<double>& std_err = path_hat.first[2];
      std_err.resize(grid_len);
      if (gain_bs.size() < 2) {
        return;
      }

      for (size_t i = 0; i < grid_len; i++) {
        // Use Welford's algorithm to get numerically stable variance estimates in one pass.
        double Mprev;
        double M;
        double Sprev = -1;
        double S;
        double n = 0;
        for (size_t b = 0; b < gain_bs.size(); b++) {
          if (gain_bs[b].size() < 1) {
            continue;
          }
          double val = gain_bs[b][i];
          if (std::isnan(val)) {
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
          // these are early grid points where min(\hat path spend) < min(path bs spend).
          std_err[i] = 0; // define these to be zero.
        }
      }

    }

    void split_sequence(std::vector<unsigned int>& result,
                        unsigned int start,
                        unsigned int end,
                        unsigned int num_parts) {
      result.reserve(num_parts + 1);

      // Return range if only 1 part
      if (num_parts == 1) {
        result.push_back(start);
        result.push_back(end + 1);
        return;
      }

      // Return vector from start to end+1 if more parts than elements
      if (num_parts > end - start + 1) {
        for (unsigned int i = start; i <= end + 1; ++i) {
          result.push_back(i);
        }
        return;
      }

      unsigned int length = (end - start + 1);
      unsigned int part_length_short = length / num_parts;
      unsigned int part_length_long = (unsigned int) std::ceil(length / ((double) num_parts));
      unsigned int cut_pos = length % num_parts;

      // Add long ranges
      for (unsigned int i = start; i < start + cut_pos * part_length_long; i = i + part_length_long) {
        result.push_back(i);
      }

      // Add short ranges
      for (unsigned int i = start + cut_pos * part_length_long; i <= end + 1; i = i + part_length_short) {
        result.push_back(i);
      }
    }

    const DataType& data;
    const SolverOptions& options;
};

} // namespace maq

#endif // MAQ_SOLVER_HPP
