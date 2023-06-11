from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.vector cimport vector

ctypedef pair[vector[vector[double]], vector[vector[size_t]]] solution_path

cdef extern from 'runner.hpp':
    solution_path run(const double* data_reward,
                      const double* data_reward_scores,
                      const double* data_cost,
                      size_t num_rows,
                      size_t num_cols,
                      double budget,
                      bool target_with_covariates,
                      unsigned int num_bootstrap,
                      unsigned int num_threads,
                      unsigned int random_seed)
