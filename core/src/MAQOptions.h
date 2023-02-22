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

#ifndef MAQ_OPTIONS_H
#define MAQ_OPTIONS_H

#include <vector>
#include <thread>

#include "Data.h"
#include "sampling/SamplingOptions.h"
#include "sampling/RandomSampler.h"

namespace maq {

typedef unsigned int uint;

struct MAQOptions {
  MAQOptions(double budget,
             size_t num_bootstrap,
             const std::vector<size_t>& clusters,
             uint samples_per_cluster,
             uint num_threads,
             uint random_seed) :
      budget(budget),
      num_bootstrap(num_bootstrap),
      clusters(clusters),
      samples_per_cluster(samples_per_cluster),
      random_seed(random_seed) {
    if (num_threads == 0) {
      num_threads = std::thread::hardware_concurrency();
    }
    this->num_threads = num_threads;
  }

  double budget;
  size_t num_bootstrap;
  const std::vector<size_t> clusters;
  uint samples_per_cluster;
  uint num_threads;
  uint random_seed;
};

} // namespace maq

#endif // MAQ_OPTIONS_H
