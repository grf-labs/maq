import cython
import numpy as np
cimport numpy as np

from libcpp cimport bool
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

from maq.maqdefs cimport solution_path, run

cpdef solver_cpp(np.ndarray[double, ndim=2, mode="c"] reward,
                 np.ndarray[double, ndim=2, mode="c"] reward_scores,
                 np.ndarray[double, ndim=2, mode="c"] cost,
                 double budget,
                 int target_with_covariates,
                 unsigned int n_bootstrap,
                 unsigned int num_threads,
                 unsigned int seed):
    cdef size_t num_rows = np.PyArray_DIMS(reward)[0]
    cdef size_t num_cols = np.PyArray_DIMS(reward)[1]

    cdef solution_path ret = run(
        &reward[0, 0], &reward_scores[0, 0], &cost[0, 0], num_rows, num_cols,
        budget, target_with_covariates, False, n_bootstrap, num_threads, seed
    )

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
        res["complete_path"] = True
    else:
        res["complete_path"] = False

    return res
