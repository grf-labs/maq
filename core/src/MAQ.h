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

#ifndef MAQ_H
#define MAQ_H

#include <vector>

#include "compute_path.h"
#include "Data.h"
#include "MAQOptions.h"
#include "sampling/SamplingOptions.h"

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
  grf::SamplingOptions sampling_options;
};

} // namespace maq

#endif // MAQ_H
