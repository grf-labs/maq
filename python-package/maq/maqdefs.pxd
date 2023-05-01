from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.vector cimport vector

ctypedef pair[vector[vector[double]], vector[vector[size_t]]] solution_path

cdef extern from 'Data.h' namespace 'maq':
    cdef cppclass Data:
        Data(const double* data_reward,
             const double* data_reward_scores,
             const double* data_cost,
             const double* data_weight,
             const int* data_tie_breaker,
             const int* clusters,
             size_t num_rows,
             size_t num_cols,
             bool col_major)

cdef extern from 'MAQOptions.h' namespace 'maq':
    cdef cppclass MAQOptions:
        MAQOptions(double budget,
                   bool paired_inference,
                   unsigned int num_bootstrap,
                   unsigned int num_threads,
                   unsigned int random_seed)

cdef extern from 'MAQ.h' namespace 'maq':
    cdef cppclass MAQ:
        MAQ(const Data& data,
            const MAQOptions& options)
        pair[solution_path, vector[vector[double]]] fit()

cdef extern from 'Data.cpp' namespace 'maq':
    pass

cdef extern from 'MAQ.cpp' namespace 'maq':
    pass

cdef extern from 'convex_hull.cpp' namespace 'maq':
    pass

cdef extern from 'compute_path.cpp' namespace 'maq':
    pass

cdef extern from 'Sampler.cpp' namespace 'maq':
    pass
