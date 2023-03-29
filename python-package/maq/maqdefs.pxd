from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.vector cimport vector

ctypedef pair[vector[vector[double]], vector[vector[size_t]]] solution_path

cdef extern from '../../core/src/Data.h' namespace 'maq':
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

cdef extern from '../../core/src/MAQOptions.h' namespace 'maq':
    cdef cppclass MAQOptions:
        MAQOptions(double budget,
                   size_t num_bootstrap,
                   unsigned int num_threads,
                   unsigned int random_seed)

cdef extern from '../../core/src/MAQ.h' namespace 'maq':
    cdef cppclass MAQ:
        MAQ(const Data& data,
            const MAQOptions& options)
        solution_path fit()

cdef extern from '../../core/src/MAQ.cpp' namespace 'maq':
    pass

cdef extern from '../../core/src/convex_hull.cpp' namespace 'maq':
    pass

cdef extern from '../../core/src/compute_path.cpp' namespace 'maq':
    pass

cdef extern from '../../core/src/Sampler.cpp' namespace 'maq':
    pass

