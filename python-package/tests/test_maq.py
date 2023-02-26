import pytest
import numpy as np
import numpy.testing as nt

from maq import MAQ
from maq.ext import solver_cpp

def test_bindings():
    budget = 100
    n = 1000
    K = 10
    reward = 1 + np.random.randn(n, K)
    cost = 0.05 + np.random.rand(n, K)
    n_bootstrap = 5
    ret = solver_cpp(reward, cost, budget, n_bootstrap, 0, 0)
    ret2 = solver_cpp(reward, cost, budget, n_bootstrap, 0, 0)

    nt.assert_equal(ret, ret2)
