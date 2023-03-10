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

#ifndef MAQ_CONVEX_HULL_H
#define MAQ_CONVEX_HULL_H

#include "Data.h"

namespace maq {

std::vector<std::vector<size_t>> convex_hull(const Data& data);

} // namespace maq

#endif // MAQ_CONVEX_HULL_H
