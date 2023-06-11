#ifndef MAQ_SAMPLER_H
#define MAQ_SAMPLER_H

#include <vector>
#include <numeric>

#include "random/random.hpp"
#include "random/algorithm.hpp"

namespace maq {

template <class DataType>
class Sampler {
  public:
  static std::vector<size_t> sample(const DataType& data,
                                    double sample_fraction,
                                    unsigned int seed) {
    std::mt19937_64 random_number_generator(seed);
    std::vector<size_t> samples;

    if (data.samples_by_cluster.empty()) {
      size_t subsample_size = static_cast<size_t>(data.get_num_rows() * sample_fraction);
      samples.resize(data.get_num_rows());
      std::iota(samples.begin(), samples.end(), 0);
      nonstd::shuffle(samples.begin(), samples.end(), random_number_generator);
      samples.resize(subsample_size);
    } else {
      // draw clusters at random
      size_t num_clusters = data.samples_by_cluster.size();
      size_t subsample_size = static_cast<size_t>(num_clusters * sample_fraction);
      std::vector<size_t> random_clusters(num_clusters);
      std::iota(random_clusters.begin(), random_clusters.end(), 0);
      nonstd::shuffle(random_clusters.begin(), random_clusters.end(), random_number_generator);
      random_clusters.resize(subsample_size);

      // fill samples-vector with samples belonging to each drawn cluster.
      for (auto cluster : random_clusters) {
        samples.insert(
          samples.end(), data.samples_by_cluster[cluster].begin(), data.samples_by_cluster[cluster].end()
        );
      }
    }

    return samples;
  }

};

} // namespace maq

#endif // MAQ_SAMPLER_H
