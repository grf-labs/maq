import numpy as np

from maq.ext import solver_cpp


class MAQ:
    """Fit a Multi-Action QINI.

    Parameters
    ----------
    n_bootstrap : int, default=200
        Number of bootstrap replicates for SEs. Default is 200.

    n_threads : int, default=0
        Number of threads used in bootstrap replicates. Default is the maximum hardware concurrency.

    seed : int, default=42
        The seed of the C++ random number generator. Default is 42.

    Attributes
    ----------
    path_spend_ : ndarray
        Fit spend path.
    path_gain_ : ndarray
        Fit gain path.
    path_std_err_ : ndarray
        Fit gain path std.err.
    path_allocated_unit_ : ndarray
        The unit allocated in spend path.
    path_allocated_arm_ : ndarray
        The arm allocated in spend path.

    Examples
    --------
    >>> import numpy as np
    >>> from maq import MAQ

    Fit a MAQ up to the maximum spend per unit.

    >>> np.random.seed(42)
    >>> n = 1000
    >>> K = 5
    >>> reward = np.random.randn(n, K)
    >>> cost = np.random.rand(n, K)
    >>> reward_eval = np.random.randn(n, K)
    >>> max_budget = np.mean(cost)

    >>> mq = MAQ()
    >>> mq.fit(reward, cost, max_budget, reward_eval)
    MAQ object.

    Get an estimate of optimal gain at a given spend along with standard errors.

    >>> mq.average_gain(spend=0.1)
    (0.005729002695991717, 0.021032360332699597)

    Get the optimal treatment allocation matrix at a given spend, a n x K array.

    >>> mq.predict(spend=0.1)
    array([[0., 0., 0., 1., 0.],
           [0., 0., 0., 0., 1.],
           [0., 0., 0., 0., 0.],
           ...,
           [0., 0., 0., 0., 1.],
           [0., 0., 0., 1., 0.],
           [0., 0., 1., 0., 0.]])

    Plot the gain curve.

    >>> import matplotlib.pyplot as plt # doctest: +SKIP

    >>> plt.plot(mq.path_spend_, mq.path_gain_) # doctest: +SKIP
    >>> plt.xlabel("Spend/unit") # doctest: +SKIP
    >>> plt.title("Gain/unit") # doctest: +SKIP
    >>> plt.show() # doctest: +SKIP
    """

    def __init__(self, n_bootstrap=200, n_threads=0, seed=42):
        assert n_threads >= 0, "n_threads should be >=0."
        assert n_bootstrap >= 0, "n_bootstrap should be >=0."
        self.n_bootstrap = n_bootstrap
        self.n_threads = n_threads
        self.seed = seed
        self._is_fit = False

    def fit(self, reward, cost, budget, reward_scores):
        """Fit the MAQ curve up to a maximum spend/user.

        Parameters
        ----------
        reward : ndarray
            A matrix of reward estimates.

        costs : ndarray
            A matrix of cost estimates.

        budget : scalar
            The maximum spend/unit to fit the MAQ path on.

        reward_scores : ndarray
            A matrix of rewards to evaluate the MAQ on.
        """

        reward = np.atleast_2d(reward)
        cost = np.atleast_2d(cost)
        reward_scores = np.atleast_2d(reward_scores)
        assert reward.shape == cost.shape, "reward and cost should have equal dims."
        assert reward.shape == reward_scores.shape, "reward and reward scores should have equal dims."
        assert np.isscalar(budget), "budget should be a scalar."
        assert (cost > 0).all(), "cost should be > 0."
        assert not np.isnan(reward).any(), "reward contains nans."
        assert not np.isnan(reward_scores).any(), "reward scores contains nans."
        assert not np.isnan(cost).any(), "cost contains nans."
        self.budget = budget

        self._path = solver_cpp(
            reward, reward_scores, cost, budget, self.n_bootstrap, self.n_threads, self.seed
        )

        self._is_fit = True
        self._dim = reward.shape
        return self

    def predict(self, spend):
        """Predict the optimal treatment allocation matrix.

        Parameters
        ----------
        spend : scalar
            The budget constraint level to predict at.

        Returns
        -------
        pi_mat : ndarray
            The optimal treatment allocation.
        """

        assert np.isscalar(spend), "spend should be a scalar."
        assert self._is_fit, "MAQ object is not fit."
        if not self._path["complete_path"]:
            assert spend <= self.budget, "maq path is not fit beyond given spend level."

        spend_grid = self._path["spend"]
        path_idx = np.searchsorted(spend_grid, spend, side="right") - 1
        pi_mat = np.zeros(self._dim, dtype="double")
        if path_idx < 0:
            return pi_mat

        ipath = self._path["ipath"][:path_idx + 1]
        kpath = self._path["kpath"][:path_idx + 1]
        ix = np.unique(ipath[::-1], return_index=True)[1]
        pi_mat[ipath[::-1][ix], kpath[::-1][ix]] = 1

        if path_idx == spend_grid.shape[0] - 1:
            return pi_mat

        # fractional adjustment?
        spend_diff = spend - spend_grid[path_idx]
        next_unit = self._path["ipath"][path_idx + 1]
        next_arm = self._path["kpath"][path_idx + 1]
        prev_arm = np.nonzero(pi_mat[next_unit, ])[0] # already assigned?

        fraction = spend_diff / (spend_grid[path_idx + 1] - spend_grid[path_idx])
        pi_mat[next_unit, next_arm] = fraction
        if prev_arm.shape[0] > 0:
            pi_mat[next_unit, prev_arm[0]] = 1 - fraction

        return pi_mat

    def average_gain(self, spend):
        """Get estimate of gain given a spend level.

        Parameters
        ----------
        spend : scalar
            The budget constraint level to predict at.

        Returns
        -------
        estimate, std_error : tuple
            Estimate of gain along with standard error.
        """

        assert np.isscalar(spend), "spend should be a scalar."
        assert self._is_fit, "MAQ object is not fit."
        if not self._path["complete_path"]:
            assert spend <= self.budget, "maq path is not fit beyond given spend level."

        spend_grid = self._path["spend"]
        path_idx = np.searchsorted(spend_grid, spend, side="right") - 1

        gain_path = self._path["gain"]
        se_path = self._path["std_err"]
        if path_idx < 0:
            estimate = 0
            std_err = 0
        elif path_idx == spend_grid.shape[0] - 1:
            estimate = gain_path[path_idx]
            std_err = se_path[path_idx]
        else:
            interp_ratio = (spend - spend_grid[path_idx]) / (spend_grid[path_idx + 1] - spend_grid[path_idx])
            estimate = gain_path[path_idx] + (gain_path[path_idx + 1] - gain_path[path_idx]) * interp_ratio
            std_err = se_path[path_idx] + (se_path[path_idx + 1] - se_path[path_idx]) * interp_ratio

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

    @property
    def path_allocated_unit_(self):
        assert self._is_fit, "MAQ object is not fit."
        return self._path["ipath"]

    @property
    def path_allocated_arm_(self):
        assert self._is_fit, "MAQ object is not fit."
        return self._path["kpath"]

    def __repr__(self):
        return "MAQ object."
