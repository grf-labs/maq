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

    mq = MAQ(n_bootstrap=n_bootstrap, n_threads=0, seed=0)
    mq.fit(reward, cost, budget)
    nt.assert_equal(ret, mq._path)

def test_MAQ():
    budget = 100
    n = 1000
    K = 10
    reward = 1 + np.random.randn(n, K)
    cost = 0.05 + np.random.rand(n, K)
    n_bootstrap = 50
    mq = MAQ(n_bootstrap=n_bootstrap, n_threads=0, seed=0)

    mq.fit(reward, cost, budget)

    nt.assert_equal(mq.average_gain(0), (0, 0))
    nt.assert_equal(mq.average_gain(100), (mq._path["gain"][-1], mq._path["std_err"][-1]))
    nt.assert_equal(mq.average_gain(mq._path["spend"][0]), (mq._path["gain"][0], mq._path["std_err"][0]))
    nt.assert_equal(mq.average_gain(mq._path["spend"][3]), (mq._path["gain"][3], mq._path["std_err"][3]))
