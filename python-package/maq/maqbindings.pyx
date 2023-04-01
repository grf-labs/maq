import cython
import numpy as np
cimport numpy as np

from libcpp cimport bool
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

from maq.maqdefs cimport Data, MAQOptions, MAQ, solution_path

cpdef solver_cpp(np.ndarray[double, ndim=2, mode="c"] reward,
                 np.ndarray[double, ndim=2, mode="c"] reward_scores,
                 np.ndarray[double, ndim=2, mode="c"] cost,
                 double budget,
                 size_t n_bootstrap,
                 unsigned int num_threads,
                 unsigned int seed):
    cdef size_t num_rows = np.PyArray_DIMS(reward)[0]
    cdef size_t num_cols = np.PyArray_DIMS(reward)[1]

    # Ignore these options for now (TODO)
    cdef double* weights_ptr = NULL
    cdef int* tie_breaker_ptr = NULL
    cdef int* clusters_ptr = NULL

    cdef bool col_major = False
    cdef Data* data_ptr
    data_ptr = new Data(
        &reward[0, 0], &reward_scores[0, 0], &cost[0, 0], weights_ptr, tie_breaker_ptr, clusters_ptr,
        num_rows, num_cols, col_major
    )

    cdef MAQOptions* opt_ptr
    opt_ptr = new MAQOptions(
        budget, n_bootstrap, num_threads, seed
    )

    cdef MAQ* maq_ptr
    maq_ptr = new MAQ(deref(data_ptr), deref(opt_ptr))
    cdef solution_path ret = deref(maq_ptr).fit()

    res = dict()
    path_len = ret.first[0].size()
    spend = np.empty(path_len, dtype="double")
    gain = np.empty(path_len, dtype="double")
    std_err = np.empty(path_len, dtype="double")
    ipath = np.empty(path_len, dtype="intp")
    kpath = np.empty(path_len, dtype="intp")

    for i in range(ret.first[0].size()):
        spend[i] = ret.first[0][i]
        gain[i] = ret.first[1][i]
        std_err[i] = ret.first[2][i]
        ipath[i] = ret.second[0][i]
        kpath[i] = ret.second[1][i]

    res["spend"] = spend
    res["gain"] = gain
    res["std_err"] = std_err
    res["ipath"] = ipath
    res["kpath"] = kpath
    if ret.second[2][0] > 0:
        res["complete_path"] = False
    else:
        res["complete_path"] = True

    del maq_ptr
    del opt_ptr
    del data_ptr

    return res
