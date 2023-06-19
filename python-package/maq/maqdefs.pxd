from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.vector cimport vector

ctypedef pair[vector[vector[double]], vector[vector[size_t]]] solution_path

cdef extern from 'wrapper.hpp':
    pair[solution_path, vector[vector[double]]] fit(
        const double* data_reward,
        const double* data_reward_scores,
        const double* data_cost,
        size_t num_rows,
        size_t num_cols,
        bool cost_matrix,
        const double* data_weight,
        const int* data_tie_breaker,
        const int* clusters,
        double budget,
        bool target_with_covariates,
        bool paired_inference,
        unsigned int num_bootstrap,
        unsigned int num_threads,
        unsigned int random_seed)
