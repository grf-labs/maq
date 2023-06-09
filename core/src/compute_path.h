#ifndef MAQ_COMPUTE_PATH_H
#define MAQ_COMPUTE_PATH_H

#include <vector>
#include "Data.h"
#include "HullData.h"

namespace maq {

typedef std::pair<std::vector<std::vector<double>>, std::vector<std::vector<size_t>>> solution_path;

solution_path compute_path(const std::vector<size_t>& samples,
                           const std::vector<std::vector<size_t>>& R,
                           const Data& data,
                           double budget,
                           bool bootstrap);

solution_path compute_path(const std::vector<size_t>& samples,
                           const std::vector<size_t>& R,
                           const HullData& data,
                           double budget,
                           bool bootstrap);

} // namespace maq

#endif // MAQ_COMPUTE_PATH_H
