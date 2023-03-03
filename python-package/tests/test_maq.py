import pytest
import doctest
import numpy as np
import numpy.testing as nt

import maq
from maq import MAQ
from maq.ext import solver_cpp

def test_bindings():
    budget = 100
    n = 1000
    K = 10
    reward = 1 + np.random.randn(n, K)
    cost = 0.05 + np.random.rand(n, K)
    n_bootstrap = 5
    ret = solver_cpp(reward, reward, cost, budget, n_bootstrap, 0, 0)
    ret2 = solver_cpp(reward, reward, cost, budget, n_bootstrap, 0, 0)

    nt.assert_equal(ret, ret2)

    mq = MAQ(n_bootstrap=n_bootstrap, n_threads=0, seed=0)
    mq.fit(reward, cost, budget, reward)
    nt.assert_equal(ret, mq._path)

def test_docstring():
    test = doctest.testmod(maq.maq)
    assert test.failed == 0

def test_MAQ():
    budget = 100
    n = 10000
    K = 10
    reward = 1 + np.random.randn(n, K)
    cost = 0.05 + np.random.rand(n, K)
    n_bootstrap = 50
    mq = MAQ(n_bootstrap=n_bootstrap, n_threads=0, seed=0)

    mq.fit(reward, cost, budget, reward)

    nt.assert_equal(mq.average_gain(0), (0, 0))
    nt.assert_equal(mq.average_gain(100), (mq._path["gain"][-1], mq._path["std_err"][-1]))
    nt.assert_equal(mq.average_gain(mq._path["spend"][0]), (mq._path["gain"][0], mq._path["std_err"][0]))
    nt.assert_equal(mq.average_gain(mq._path["spend"][3]), (mq._path["gain"][3], mq._path["std_err"][3]))

    nt.assert_almost_equal(
        mq.average_gain(100)[0],
        np.sum(mq.predict(100) * reward) / n,
        decimal=10
    )
    nt.assert_almost_equal(
        mq.average_gain(0)[0],
        np.sum(mq.predict(0) * reward) / n,
        decimal=10
    )
    nt.assert_almost_equal(
        mq.average_gain(0.1)[0],
        np.sum(mq.predict(0.1) * reward) / n,
        decimal=10
    )
    nt.assert_almost_equal(
        mq.average_gain(0.25)[0],
        np.sum(mq.predict(0.25) * reward) / n,
        decimal=10
    )
    sp = mq._path["gain"][0]
    nt.assert_almost_equal(
        mq.average_gain(sp)[0],
        np.sum(mq.predict(sp) * reward) / n,
        decimal=10
    )
    sp = mq._path["gain"][5]
    nt.assert_almost_equal(
        mq.average_gain(sp)[0],
        np.sum(mq.predict(sp) * reward) / n,
        decimal=10
    )
    sp = mq._path["gain"][-1]
    nt.assert_almost_equal(
        mq.average_gain(sp)[0],
        np.sum(mq.predict(sp) * reward) / n,
        decimal=10
    )
    sp = mq._path["gain"][-2]
    nt.assert_almost_equal(
        mq.average_gain(sp)[0],
        np.sum(mq.predict(sp) * reward) / n,
        decimal=10
    )

    reward2 = np.random.randn(n, 1)
    cost2 = np.random.rand(n, 1)
    mq.fit(reward2, cost2, budget, reward2)
    nt.assert_equal(
        np.nonzero(mq.predict(budget)),
        np.where(reward2 > 0)
    )
