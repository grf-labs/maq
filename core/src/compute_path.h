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

#ifndef MAQ_COMPUTE_PATH_H
#define MAQ_COMPUTE_PATH_H

#include <vector>
#include "Data.h"

namespace maq {

typedef std::pair<std::vector<std::vector<double>>, std::vector<std::vector<size_t>>> solution_path;

bool equal_doubles(double first, double second, double epsilon);

solution_path compute_path(const std::vector<size_t>& samples,
                           const std::vector<std::vector<size_t>>& R,
                           const Data& data,
                           double budget,
                           bool bootstrap);
} // namespace maq

#endif // MAQ_COMPUTE_PATH_H
