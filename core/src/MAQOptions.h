#ifndef MAQ_OPTIONS_H
#define MAQ_OPTIONS_H

#include <vector>
#include <thread>

#include "Data.h"

namespace maq {

struct MAQOptions {
  MAQOptions(double budget,
             unsigned int num_bootstrap,
             unsigned int num_threads,
             unsigned int random_seed) :
      budget(budget),
      num_bootstrap(num_bootstrap),
      random_seed(random_seed) {
    if (num_threads == 0) {
      num_threads = std::thread::hardware_concurrency();
    }
    this->num_threads = num_threads;
  }

  double budget;
  unsigned int num_bootstrap;
  unsigned int num_threads;
  unsigned int random_seed;
};

} // namespace maq

#endif // MAQ_OPTIONS_H
