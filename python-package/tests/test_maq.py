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
    ret = solver_cpp(reward, reward, cost, budget, True, n_bootstrap, True, 0, 0)
    ret2 = solver_cpp(reward, reward, cost, budget, True, n_bootstrap, True, 0, 0)
    nt.assert_equal(ret, ret2)

    mq = MAQ(budget, n_bootstrap=n_bootstrap, n_threads=0, seed=0)
    mq.fit(reward, cost, reward)
    nt.assert_equal(ret, mq._path)

    mq_avg = MAQ(
        budget,
        target_with_covariates=False,
        n_bootstrap=n_bootstrap,
        n_threads=0,
        seed=0,
    )
    mq_avg.fit(reward, cost, reward)
    nt.assert_almost_equal(
        mq_avg.average_gain(100)[0],
        np.sum(mq_avg.predict(100) * reward) / n,
        decimal=10,
    )


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
    mq = MAQ(budget, n_bootstrap=n_bootstrap, n_threads=0, seed=0)

    mq.fit(reward, cost, reward)

    nt.assert_equal(mq.average_gain(0), (0, 0))
    nt.assert_equal(
        mq.average_gain(100), (mq._path["gain"][-1], mq._path["std_err"][-1])
    )
    nt.assert_equal(
        mq.average_gain(mq._path["spend"][0]),
        (mq._path["gain"][0], mq._path["std_err"][0]),
    )
    nt.assert_equal(
        mq.average_gain(mq._path["spend"][3]),
        (mq._path["gain"][3], mq._path["std_err"][3]),
    )

    nt.assert_almost_equal(
        mq.average_gain(100)[0], np.sum(mq.predict(100) * reward) / n, decimal=10
    )
    nt.assert_almost_equal(
        mq.average_gain(0)[0], np.sum(mq.predict(0) * reward) / n, decimal=10
    )
    nt.assert_almost_equal(
        mq.average_gain(0.1)[0], np.sum(mq.predict(0.1) * reward) / n, decimal=10
    )
    nt.assert_almost_equal(
        mq.average_gain(0.25)[0], np.sum(mq.predict(0.25) * reward) / n, decimal=10
    )
    sp = mq._path["gain"][0]
    nt.assert_almost_equal(
        mq.average_gain(sp)[0], np.sum(mq.predict(sp) * reward) / n, decimal=10
    )
    sp = mq._path["gain"][5]
    nt.assert_almost_equal(
        mq.average_gain(sp)[0], np.sum(mq.predict(sp) * reward) / n, decimal=10
    )
    sp = mq._path["gain"][-1]
    nt.assert_almost_equal(
        mq.average_gain(sp)[0], np.sum(mq.predict(sp) * reward) / n, decimal=10
    )
    sp = mq._path["gain"][-2]
    nt.assert_almost_equal(
        mq.average_gain(sp)[0], np.sum(mq.predict(sp) * reward) / n, decimal=10
    )

    reward2 = np.random.randn(n, 1)
    cost2 = np.random.rand(n, 1)
    mq.fit(reward2, cost2, reward2)
    nt.assert_equal(np.nonzero(mq.predict(budget)), np.where(reward2 > 0))


def test_cost_arg():
    n = 1000
    K = 2
    reward = 1 + np.random.randn(n, K)
    cost_vec = np.array([1, 2])
    cost = np.repeat(cost_vec[None, :], n, axis=0)

    mq = MAQ().fit(reward, cost, reward)
    mq2 = MAQ().fit(reward, cost_vec, reward)
    mq3 = MAQ().fit(reward, [1, 2], reward)

    nt.assert_equal(mq._path, mq2._path)
    nt.assert_equal(mq._path, mq3._path)

    mq4 = MAQ().fit(reward[:, 0], cost[:, 0], reward[:, 0])
    mq5 = MAQ().fit(reward[:, 0], cost_vec[0], reward[:, 0])
    mq6 = MAQ().fit(reward[:, 0], 1, reward[:, 0])

    nt.assert_equal(mq4._path, mq5._path)
    nt.assert_equal(mq4._path, mq6._path)


def test_prediction_type():
    n = 1000
    K = 10
    reward = 1 + np.random.randn(n, K)
    cost = 0.05 + np.random.rand(n, K)

    mq = MAQ().fit(reward, cost, reward)
    sp = mq.path_spend_[50]
    pi_mat = mq.predict(sp)
    nt.assert_almost_equal(np.sum(pi_mat * cost) / n, sp, decimal=10)
    pi_vec = mq.predict(sp, prediction_type="vector")
    csum = 0
    for i, k in enumerate(pi_vec):
        if k == 0:
            continue
        csum = csum + cost[i, k - 1]
    nt.assert_almost_equal(csum / n, sp, decimal=10)


def test_difference_gain():
    # these gains should be statistically indistinguishable with 95 % coverage
    res = []
    spend = 0.2
    for _ in range(250):
        n = 1000
        K = 3
        reward = np.random.rand(n, K)
        cost = np.random.rand(n, K)
        reward_eval = np.random.randn(n, K)
        mq1 = MAQ(n_bootstrap=200).fit(reward, cost, reward_eval)
        mq2 = MAQ(n_bootstrap=200).fit(
            reward + np.random.randn(n)[:, None], cost, reward_eval
        )
        est, sd = mq1.difference_gain(mq2, spend)

        coverage = int(abs(est / sd) <= 1.96)
        res.append(coverage)

    nt.assert_allclose(np.mean(res), 0.95, rtol=0.05)

    res_avg = []
    for _ in range(250):
        n = 1000
        K = 3
        reward = np.random.rand(n, K)
        cost = np.random.rand(n, K)
        reward_eval = np.random.randn(n, K)

        mq1 = MAQ(n_bootstrap=200).fit(reward, cost, reward_eval)
        mq2 = MAQ(n_bootstrap=200, target_with_covariates=False).fit(
            reward, cost, reward_eval
        )
        est, sd = mq1.difference_gain(mq2, spend)

        coverage = int(abs(est / sd) <= 1.96)
        res_avg.append(coverage)

    nt.assert_allclose(np.mean(res_avg), 0.95, rtol=0.05)


def test_integrated_difference_gain():
    # these gains should be statistically indistinguishable with 95 % coverage
    res = []
    spend = 0.2
    for _ in range(250):
        n = 1000
        K = 3
        reward = np.random.rand(n, K)
        cost = np.random.rand(n, K)
        reward_eval = np.random.randn(n, K)
        mq1 = MAQ(n_bootstrap=200).fit(reward, cost, reward_eval)
        mq2 = MAQ(n_bootstrap=200).fit(
            reward + np.random.randn(n)[:, None], cost, reward_eval
        )
        est, sd = mq1.integrated_difference(mq2, spend)

        coverage = int(abs(est / sd) <= 1.96)
        res.append(coverage)

    nt.assert_allclose(np.mean(res), 0.95, rtol=0.05)

    res_avg = []
    for _ in range(250):
        n = 1000
        K = 3
        reward = np.random.rand(n, K)
        cost = np.random.rand(n, K)
        reward_eval = np.random.randn(n, K)

        mq1 = MAQ(n_bootstrap=200).fit(reward, cost, reward_eval)
        mq2 = MAQ(n_bootstrap=200, target_with_covariates=False).fit(
            reward, cost, reward_eval
        )
        est, sd = mq1.integrated_difference(mq2, spend)

        coverage = int(abs(est / sd) <= 1.96)
        res_avg.append(coverage)

    nt.assert_allclose(np.mean(res_avg), 0.95, rtol=0.05)


def test_integrated_difference_grid_numerics():
    n = 10
    K = 1
    reward1 = np.random.rand(n, K)
    reward2 = np.random.rand(n, K)
    cost = 1
    reward_eval = np.random.randn(n, K)
    q1 = MAQ().fit(reward1, cost, reward_eval)
    q2 = MAQ().fit(reward2, cost, reward_eval)
    sp = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    g1 = 0
    g2 = 0
    for i, spend in enumerate(sp):
        g1 += q1.average_gain(spend)[0]
        g2 += q2.average_gain(spend)[0]

        est1 = (g1 - g2) / (i + 1)
        est2 = np.mean(q1.path_gain_[: i + 1]) - np.mean(q2.path_gain_[: i + 1])
        est = q1.integrated_difference(q2, spend)[0]
        nt.assert_allclose(est1, est2, rtol=1e-07)
        nt.assert_allclose(est, est2, rtol=1e-07)


def test_basic_perf_timing():
    import time

    n = 1000000
    K = 5
    reward = np.random.rand(n, K)
    cost = np.random.rand(n, K)
    reward_eval = np.random.randn(n, K)

    mq = MAQ()
    start = time.time()
    [mq.fit(reward, cost, reward_eval) for _ in range(4)]
    end = time.time()
    elapsed = end - start
    nt.assert_array_less(elapsed / 4, 2.7)
