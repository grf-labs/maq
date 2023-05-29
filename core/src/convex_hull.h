#ifndef MAQ_CONVEX_HULL_H
#define MAQ_CONVEX_HULL_H

#include "HullData.h"

namespace maq {

std::vector<std::vector<size_t>> convex_hull(const HullData& data);

} // namespace maq

#endif // MAQ_CONVEX_HULL_H
