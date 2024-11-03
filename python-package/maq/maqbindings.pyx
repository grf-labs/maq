import cython
import numpy as np
cimport numpy as np
from libcpp cimport bool

from maq.maqdefs cimport pair, vector, solution_path, run

cpdef solver_cpp(np.ndarray[double, ndim=2, mode="c"] reward,
                 np.ndarray[double, ndim=2, mode="c"] reward_scores,
                 np.ndarray[double, ndim=2, mode="c"] cost,
                 double budget,
                 int target_with_covariates,
                 unsigned int n_bootstrap,
                 int paired_inference,
                 unsigned int num_threads,
                 unsigned int seed):
    cdef size_t num_rows = np.PyArray_DIMS(reward)[0]
    cdef size_t num_cols = np.PyArray_DIMS(reward)[1]
    cdef bool cost_matrix = True
    if cost.shape[0] == 1:
        cost_matrix = False

    # Ignore these options for now (TODO)
    cdef double* weights_ptr = NULL
    cdef int* tie_breaker_ptr = NULL
    cdef int* clusters_ptr = NULL

    cdef pair[solution_path, vector[vector[double]]] ret = run(
        &reward[0, 0],
        &reward_scores[0, 0],
        &cost[0, 0],
        num_rows,
        num_cols,
        cost_matrix,
        weights_ptr,
        tie_breaker_ptr,
        clusters_ptr,
        budget,
        target_with_covariates,
        paired_inference,
        n_bootstrap,
        num_threads,
        seed
    )
    cdef solution_path path = ret.first

    res = dict()
    path_len = path.first[0].size()
    spend = np.empty(path_len, dtype="double")
    gain = np.empty(path_len, dtype="double")
    std_err = np.empty(path_len, dtype="double")
    ipath = np.empty(path_len, dtype="int")
    kpath = np.empty(path_len, dtype="int")
    if paired_inference:
        gain_bs = np.empty((n_bootstrap, path_len), dtype="double")
    else:
        gain_bs = np.empty((0, 0), dtype="double")
    # faster copy into nparrays with memoryviews
    cdef double[::] view_spend = spend
    cdef double[::] view_gain = gain
    cdef double[::] view_std_err = std_err
    cdef int[::] view_ipath = ipath
    cdef int[::] view_kpath = kpath
    cdef double[:, ::1] view_gain_bs = gain_bs

    for i in range(path_len):
        view_spend[i] = path.first[0][i]
        view_gain[i] = path.first[1][i]
        view_std_err[i] = path.first[2][i]
        view_ipath[i] = path.second[0][i]
        view_kpath[i] = path.second[1][i]
    if paired_inference:
        for b in range(n_bootstrap):
            for i in range(path_len):
                view_gain_bs[b, i] = ret.second[b][i]

    res["spend"] = spend
    res["gain"] = gain
    res["std_err"] = std_err
    res["ipath"] = ipath
    res["kpath"] = kpath
    if path.second[2][0] > 0:
        res["complete_path"] = True
    else:
        res["complete_path"] = False
    res["gain_bs"] = gain_bs

    return res
