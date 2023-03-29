#ifndef MAQ_SAMPLER_H
#define MAQ_SAMPLER_H

#include <vector>
#include "Data.h"

namespace maq {

class Sampler {
public:
  static std::vector<size_t> sample(const Data& data, double sample_fraction, unsigned int seed);
};

} // namespace maq

#endif // MAQ_SAMPLER_H
