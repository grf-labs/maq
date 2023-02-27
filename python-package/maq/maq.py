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
        self.budget = budget

        self._path = solver_cpp(reward, cost, budget,
            self.n_bootstrap, self.n_threads, self.seed)

        self._is_fit = True
        self._dim = reward.shape
        return self

    def predict(self, spend):
        assert np.isscalar(spend), "spend should be a scalar."
        assert self._is_fit, "MAQ object is not fit."
        if not self._path["complete_path"]:
            assert spend <= self.budget, "maq path is not fit beyond given spend level."

        spend_grid = self._path["spend"]
        path_idx = np.searchsorted(spend_grid, spend, side = "right") - 1
        if path_idx < 0:
            return np.zeros(self._dim, dtype = "intp")


        ipath = self._path["ipath"]
        kpath = self._path["kpath"]

        return 42

    def average_gain(self, spend):
        assert np.isscalar(spend), "spend should be a scalar."
        assert self._is_fit, "MAQ object is not fit."
        if not self._path["complete_path"]:
            assert spend <= self.budget, "maq path is not fit beyond given spend level."

        spend_grid = self._path["spend"]
        path_idx = np.searchsorted(spend_grid, spend, side = "right") - 1

        gain_path = self._path["gain"]
        se_path = self._path["std_err"]
        if path_idx < 0:
            estimate = 0
            std_err = 0
        elif path_idx == spend_grid.shape[0] - 1:
            estimate = gain_path[path_idx]
            std_err = se_path[path_idx]
        else:
            interp_ratio = (spend - spend_grid[path_idx]) / (spend_grid[path_idx+1] - spend_grid[path_idx])
            estimate = gain_path[path_idx] + (gain_path[path_idx+1] - gain_path[path_idx]) * interp_ratio
            std_err = se_path[path_idx] + (se_path[path_idx+1] - se_path[path_idx]) * interp_ratio

        return estimate, std_err


    @property
    def path_spend_(self):
        assert self._is_fit, "MAQ object is not fit."
        return self._path["spend"]

    @property
    def path_gain_(self):
        assert self._is_fit, "MAQ object is not fit."
        return self._path["gain"]

    @property
    def path_std_err_(self):
        assert self._is_fit, "MAQ object is not fit."
        return self._path["std_err"]

    def __repr__(self):
        return "MAQ object."
