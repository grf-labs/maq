// Multi-Action Qini (maq).
// https://github.com/grf-labs/maq
// Distributed under the MIT License.

#ifndef MAQ_H
#define MAQ_H

#include <vector>

#include "compute_path.h"
#include "Data.h"
#include "MAQOptions.h"

namespace maq {

typedef unsigned int uint;

class MAQ {
  public:
  MAQ(const Data& data,
      const MAQOptions& options);

  solution_path fit();

  private:
  std::vector<std::vector<double>> fit_paths(const solution_path& path_hat,
                                             const std::vector<std::vector<size_t>>& R);

  std::vector<std::vector<double>> fit_paths_batch(size_t start_index,
                                                   size_t num_replicates,
                                                   const solution_path& path_hat,
                                                   const std::vector<std::vector<size_t>>& R);

  std::vector<double> interpolate_path(const solution_path& path_hat,
                                       const solution_path& path_hat_b);

  void compute_std_err(solution_path& path_hat,
                       const std::vector<std::vector<double>>& gain_bs);

  void split_sequence(std::vector<uint>& result, uint start, uint end, uint num_parts);

  Data data;
  MAQOptions options;
};

} // namespace maq

#endif // MAQ_H
