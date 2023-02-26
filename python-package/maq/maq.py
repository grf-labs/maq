import numpy as np

from maq.ext import solver_cpp


class MAQ:
    r"""
    Fit a Multi-Action QINI.


    """

    def __init__(self, n_bootstrap=200, n_threads=0, seed=42):
        assert n_threads >= 0, "n_threads should be >=0."
        assert n_bootstrap >= 0, "n_bootstrap should be >=0."
        self.n_bootstrap = n_bootstrap
        self.n_threads = n_threads
        self.seed = seed
        self._is_fit = False

    def fit(self, reward, cost, budget):
        reward = np.atleast_2d(reward)
        cost = np.atleast_2d(cost)
        assert reward.shape == cost.shape, "reward and cost should have equal dims."
        assert np.isscalar(budget), "budget should be a scalar."
        assert (cost > 0).all(), "cost should be > 0."
        assert not np.isnan(reward).any(), "reward contains nans."
        assert not np.isnan(cost).any(), "cost contains nans."

        self._path = solver_cpp(reward, cost, budget,
            self.n_bootstrap, self.n_threads, self.seed)

        self._is_fit = True
        return self

    def predict(spend):

        return 1

    def average_gain(spend):

        return 1
