import numpy as np

from maq.ext import solver_cpp


def get_ipw_scores(Y, W, W_hat=None):
    """Construct evaluation scores via inverse-propensity weighting.

    Parameters
    ----------
    Y : ndarray
        A vector of test set outcomes.

    W : ndarray
        A vector of test set treatment assignments encoded as integers 0,...,K where
        k=0 is the control arm and k=1,..,K one of K treatment arms.

    W_hat : ndarray, default=None
        Optional num_samples * (K + 1) array of treatment propensities where the k-th column
        contains the propensity scores for the k-th arm. If None (Default), then the assignment
        probabilities are assumed to be uniform and the same for each arm.

    Returns
    -------
    ndarray
        An (num_samples * K) array of scores.

    Examples
    --------
    Generate some synthetic data from K=2 treatment arms with equal assignment probabilities.

    >>> import numpy as np
    >>> np.random.seed(42)
    >>> n = 2000
    >>> W = np.random.choice([0, 1, 2], n)
    >>> Y = 42 * (W == 1) - 42 * (W == 2) + np.random.rand(n)
    >>> IPW_scores = get_ipw_scores(Y, W)

    An IPW estimate of E[Y(1) - Y(0)] and E[Y(2) - Y(0)] ~ 42 and -42.

    >>> IPW_scores.mean(axis=0)
    array([ 41.51013928, -41.09905001])

    Draw non-uniformly from the different arms.

    >>> p = np.array([0.2, 0.2, 0.6])
    >>> W = np.random.choice([0, 1, 2], n, p = p)
    >>> Y = 42 * (W == 1) - 42 * (W == 2) + np.random.rand(n)
    >>> # Construct a matrix of treatment assignment probabilities.
    >>> W_hat = np.repeat(p[None,:], n, axis=0)
    >>> IPW_scores = get_ipw_scores(Y, W, W_hat)
    >>> IPW_scores.mean(axis=0)
    array([ 43.06614088, -41.60718507])
    """
    if Y.ndim > 1 and Y.shape != W.shape:
        raise ValueError("Y and W should be equal length vectors.")
    K = max(W)
    if (
        not isinstance(K, np.integer)
        or K <= 0
        or not np.array_equal(np.unique(W), np.array(range(K + 1)))
    ):
        raise ValueError("W should be a vector of integers (0, 1, ..., K).")

    if W_hat is None:
        W_hat = 1.0 / (K + 1)
    elif W_hat.shape != (len(W), K + 1):
        raise ValueError("W_hat should be a matrix of propensities.")
    elif (
        np.any(W_hat <= 0)
        or np.any(W_hat >= 1)
        or np.any(abs(W_hat.sum(axis=1) - 1) > 1e-6)
    ):
        raise ValueError("W_hat entries should be in (0, 1) and sum to 1.")

    Y_mat = np.zeros((len(W), K + 1))
    Y_mat[range(len(W)), W] = Y
    Y_ipw = Y_mat / W_hat

    return Y_ipw[:, 1:] - Y_ipw[:, 0][:, None]


class MAQ:
    """Fit a Multi-Armed Qini.

    Given n test set samples and K treatment arms, construct a Qini curve Q(B) that quantifies the
    value of assigning treatment in accordance with an estimated treatment effect function
    while satisfying a budget constraint B.

    Parameters
    ----------
    budget : scalar, default=None
        The maximum spend per unit to fit the Qini curve on.
        Setting this to None (Default), will fit the path up to a maximum spend per unit
        where each unit that is expected to benefit is treated (that is, hat tau_k(X_i) > 0).

    target_with_covariates : bool, default=True
        If TRUE, then the optimal policy takes covariates into account. If FALSE, then the optimal policy
        only takes the average reward and cost into account when allocating treatment. Can be used to
        construct a baseline Qini curve to assess the value of targeting with covariates.

    n_bootstrap : int, default=0
        Number of bootstrap replicates for SEs. Default is 0.

    paired_inference : bool, default=True
        Whether to allow for paired tests with other Qini curves fit on the same evaluation data.
        If TRUE (Default) then the path of bootstrap replicates are stored in order to perform
        paired comparisons that account for the correlation between curves evaluated on the same data.
        This takes memory on the order of O(n_bootstrap*num_samples*K) and requires the comparison
        objects to be fit with the same seed and n_bootstrap values as well as the same number of samples.

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

    Fit a Qini curve on toy data.

    >>> np.random.seed(42)
    >>> n = 1000
    >>> K = 5
    >>> tau_hat = np.random.randn(n, K)
    >>> cost = np.random.rand(n, K)
    >>> DR_scores = np.random.randn(n, K)

    >>> mq = MAQ(n_bootstrap=200)
    >>> mq.fit(tau_hat, cost, DR_scores)
    MAQ object with 1000 units and 5 arms.

    Get an estimate of gain at a given spend along with standard errors.

    >>> mq.average_gain(spend=0.1)
    (0.005729002695991717, 0.019814651108894354)

    Get the underlying treatment allocation matrix at a given spend, a n x K array.

    >>> mq.predict(spend=0.1)
    array([[0., 0., 0., 1., 0.],
           [0., 0., 0., 0., 1.],
           [0., 0., 0., 0., 0.],
           ...,
           [0., 0., 0., 0., 1.],
           [0., 0., 0., 1., 0.],
           [0., 0., 1., 0., 0.]])

    Plot the Qini curve (requires matplotlib).

    >>> mq.plot() # doctest: +SKIP

    Show 95% confidence bars.

    >>> mq.plot(show_ci=True) # doctest: +SKIP
    """

    def __init__(
        self,
        budget=None,
        target_with_covariates=True,
        n_bootstrap=0,
        paired_inference=True,
        n_threads=0,
        seed=42,
    ):
        if budget is None:
            budget = np.finfo(np.float64).max
        assert np.isscalar(budget), "budget should be a scalar."
        assert n_threads >= 0, "n_threads should be >=0."
        assert n_bootstrap >= 0, "n_bootstrap should be >=0."
        self.budget = budget
        self.target_with_covariates = target_with_covariates
        self.n_bootstrap = n_bootstrap
        self.paired_inference = paired_inference
        self.n_threads = n_threads
        self.seed = seed
        self._is_fit = False

    def fit(self, reward, cost, DR_scores):
        """Fit a Qini curve.

        Parameters
        ----------
        reward : ndarray
            A matrix of reward estimates with rows corresponding to units and columns containing
            the treatment effect estimates for the K different arms.

        cost : ndarray
            A matrix of costs. If the costs only vary by arm and not by unit, then this
            can also be a K-vector of costs for each arm.

        DR_scores : ndarray
            A matrix of evaluation scores to estimate the Qini curve on.
        """
        # ensure dims are (n, K)
        reward = np.reshape(reward, (reward.shape[0], -1))
        DR_scores = np.reshape(DR_scores, (DR_scores.shape[0], -1))
        # if costs are the same for each unit/arm, they can have dim (1, K)
        if len(np.atleast_1d(cost)) == reward.shape[1]:
            cost = np.atleast_2d(cost).astype(float)
        else:
            cost = np.reshape(cost, (np.atleast_1d(cost).shape[0], -1)).astype(float)

        if reward.shape != DR_scores.shape or cost.shape[1] != reward.shape[1]:
            raise ValueError(
                "reward, costs, and evaluation scores should have conformable dimensions."
            )
        if cost.shape[0] > 1 and cost.shape[0] != reward.shape[0]:
            raise ValueError(
                "reward, costs, and evaluation scores should have conformable dimensions."
            )
        if np.any(cost <= 0):
            raise ValueError("cost should be > 0.")

        if np.isnan(reward).any() or np.isnan(DR_scores).any() or np.isnan(cost).any():
            raise ValueError(
                "reward, costs, and evaluation scores should have no missing values."
            )

        self._path = solver_cpp(
            np.ascontiguousarray(reward),
            np.ascontiguousarray(DR_scores),
            np.ascontiguousarray(cost),
            self.budget,
            self.target_with_covariates,
            self.n_bootstrap,
            self.paired_inference,
            self.n_threads,
            self.seed,
        )

        self._is_fit = True
        self._dim = reward.shape
        return self

    def predict(self, spend, prediction_type="matrix"):
        """Predict treatment allocation.

        Parameters
        ----------
        spend : scalar
            The budget constraint level to predict at.

        type : str
            If "matrix", then represent the underlying treatment allocation as a num_samples * K
            matrix, where for row i, the k-th element is 1 if assigning the k-th arm to unit i is
            optimal at a given spend, and 0 otherwise (with all entries 0 if the control arm is assigned).
            If "vector" then represent the underlying treatment allocation as a num_samples-length
            vector where entries take values in the set k=(0, 1, ..., K) where k=0 represents the
            control arm. If at a given spend, the treatment allocation is fractional (i.e., for a
            single unit i there is not sufficient budget left to assign the i-th unit an initial arm,
            or upgrade to the next costlier arm), then the returned vector is the treatment
            allocation in the solution path where the allocation is integer-valued but incurs a cost
            (slightly) less than 'spend'.

        Returns
        -------
        pi_mat : ndarray
            The treatment allocation at a given spend per unit.
        """

        assert np.isscalar(spend), "spend should be a scalar."
        self._ensure_fit()
        if not self._path["complete_path"] and spend > self.budget:
            raise ValueError("maq path is not fit beyond given spend level.")

        spend_grid = self._path["spend"]
        path_idx = np.searchsorted(spend_grid, spend, side="right") - 1
        if path_idx < 0:
            if prediction_type == "matrix":
                return np.zeros(self._dim, dtype="double")
            else:
                return np.zeros(self._dim[0], dtype="int")

        ipath = self._path["ipath"][: path_idx + 1]
        kpath = self._path["kpath"][: path_idx + 1]
        ix = np.unique(ipath[::-1], return_index=True)[1]

        if prediction_type == "vector":
            pi_vec = np.zeros(self._dim[0], dtype="int")
            pi_vec[ipath[::-1][ix]] = kpath[::-1][ix] + 1
            return pi_vec

        pi_mat = np.zeros(self._dim, dtype="double")
        pi_mat[ipath[::-1][ix], kpath[::-1][ix]] = 1

        if path_idx == spend_grid.shape[0] - 1:
            return pi_mat

        # fractional adjustment?
        spend_diff = spend - spend_grid[path_idx]
        next_unit = self._path["ipath"][path_idx + 1]
        next_arm = self._path["kpath"][path_idx + 1]
        prev_arm = np.nonzero(pi_mat[next_unit,])[0]  # already assigned?

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
            The spend level.

        Returns
        -------
        estimate, std_error : tuple
            Estimate of gain along with standard errors.
        """

        assert np.isscalar(spend), "spend should be a scalar."
        self._ensure_fit()
        if not self._path["complete_path"] and spend > self.budget:
            raise ValueError("maq path is not fit beyond given spend level.")

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
            interp_ratio = (spend - spend_grid[path_idx]) / (
                spend_grid[path_idx + 1] - spend_grid[path_idx]
            )
            estimate = (
                gain_path[path_idx]
                + (gain_path[path_idx + 1] - gain_path[path_idx]) * interp_ratio
            )
            std_err = (
                se_path[path_idx]
                + (se_path[path_idx + 1] - se_path[path_idx]) * interp_ratio
            )

        return estimate, std_err

    def difference_gain(self, other, spend):
        """Get estimate of difference in gain given a spend level with paired standard errors.

        Parameters
        ----------
        other : MAQ object.
            The other Qini curve to subtract with.

        spend : scalar
            The spend level.

        Returns
        -------
        estimate, std_error : tuple
            Estimate of difference in gain along with standard errors.
        """
        assert np.isscalar(spend), "spend should be a scalar."
        self._ensure_fit()
        if not self._path["complete_path"] and spend > self.budget:
            raise ValueError("maq path is not fit beyond given spend level.")
        other._ensure_fit()
        if not other._path["complete_path"] and spend > other.budget:
            raise ValueError("comparison maq path is not fit beyond given spend level.")
        if (
            self.seed != other.seed
            or self.n_bootstrap != other.n_bootstrap
            or self._dim[0] != other._dim[0]
            or not self.paired_inference
            or not other.paired_inference
        ):
            raise ValueError(
                """Paired comparisons require maq objects to be fit with paired_inference=True
                as well as with the same random seed, bootstrap replicates, and data."""
            )
        point_estimate = self.average_gain(spend)[0] - other.average_gain(spend)[0]
        if self.n_bootstrap < 2:
            return point_estimate, 0

        # Compute paired std.errors
        def _get_estimates(obj):
            gain_bs = obj._path["gain_bs"]
            spend_grid = obj._path["spend"]
            path_idx = np.searchsorted(spend_grid, spend, side="right") - 1
            if path_idx < 0:
                estimates = 0
            elif path_idx == spend_grid.shape[0] - 1:
                estimates = gain_bs[:, path_idx]
            else:
                interp_ratio = (spend - spend_grid[path_idx]) / (
                    spend_grid[path_idx + 1] - spend_grid[path_idx]
                )
                estimates = (
                    gain_bs[:, path_idx]
                    + (gain_bs[:, path_idx + 1] - gain_bs[:, path_idx]) * interp_ratio
                )
            return estimates

        estimates_lhs = _get_estimates(self)
        estimates_rhs = _get_estimates(other)
        std_err = np.nanstd(estimates_lhs - estimates_rhs)
        if np.isnan(std_err):
            std_err = 0

        return point_estimate, std_err

    def integrated_difference(self, other, spend):
        """Get estimate of the area between two Qini curves with paired standard errors.

        Parameters
        ----------
        other : MAQ object.
            The other Qini curve to subtract with.

        spend : scalar
            The spend level.

        Returns
        -------
        estimate, std_error : tuple
            An estimate of the area between the two curves along with standard errors.
        """
        assert np.isscalar(spend), "spend should be a scalar."
        self._ensure_fit()
        if not self._path["complete_path"] and spend > self.budget:
            raise ValueError("maq path is not fit beyond given spend level.")
        other._ensure_fit()
        if not other._path["complete_path"] and spend > other.budget:
            raise ValueError("comparison maq path is not fit beyond given spend level.")
        if (
            self.seed != other.seed
            or self.n_bootstrap != other.n_bootstrap
            or self._dim[0] != other._dim[0]
            or not self.paired_inference
            or not other.paired_inference
        ):
            raise ValueError(
                """Paired comparisons require maq objects to be fit with paired_inference=True
                as well as with the same random seed, bootstrap replicates, and data."""
            )

        # Estimate an AUC via estimating the difference \int_{0}^{\bar B} Q_a(B)dB - \int_{0}^{\bar B} Q_b(B)dB.
        def _get_estimates(gain, spend_grid, path_idx):
            if path_idx < 0:
                estimates = np.asarray([0])
            elif path_idx == spend_grid.shape[0] - 1:
                area_offset = 0
                # Are we summing beyond the point at which the curve plateaus?
                if spend > spend_grid[path_idx]:
                    spend_delta = spend - spend_grid[path_idx]
                    area_offset = gain[:, path_idx] * spend_delta
                estimates = np.nanmean(gain, axis=1) + area_offset
            else:
                interp_ratio = (spend - spend_grid[path_idx]) / (
                    spend_grid[path_idx + 1] - spend_grid[path_idx]
                )
                if interp_ratio < 1e-10:
                    adj = np.repeat(np.nan, gain.shape[0])
                else:
                    adj = (
                        gain[:, path_idx]
                        + (gain[:, path_idx + 1] - gain[:, path_idx]) * interp_ratio
                    )
                estimates = np.nanmean(
                    np.hstack((gain[:, : path_idx + 1], adj[:, None])), axis=1
                )

            return estimates

        path_idx_lhs = np.searchsorted(self._path["spend"], spend, side="right") - 1
        path_idx_rhs = np.searchsorted(other._path["spend"], spend, side="right") - 1

        point_estimate = (
            _get_estimates(
                self._path["gain"][None, :], self._path["spend"], path_idx_lhs
            )[0]
            - _get_estimates(
                other._path["gain"][None, :], other._path["spend"], path_idx_rhs
            )[0]
        )
        if self.n_bootstrap < 2:
            return point_estimate, 0
        # Compute paired std.errors
        estimates_lhs = _get_estimates(
            self._path["gain_bs"], self._path["spend"], path_idx_lhs
        )
        estimates_rhs = _get_estimates(
            other._path["gain_bs"], other._path["spend"], path_idx_rhs
        )
        std_err = np.nanstd(estimates_lhs - estimates_rhs)
        if np.isnan(std_err):
            std_err = 0

        return point_estimate, std_err

    def plot(self, show_ci=False, **kwargs):
        """Plot the Qini curve (requires matplotlib).

        If the underlying policy involves treating zero units (as would be the case if all
        reward estimates are negative or the average is <0), then nothing is plot.

        Parameters
        ----------
        show_ci : bool
            Whether to show estimated 95% confidence bars.
        **kwargs : additional arguments passed to matplotlib.pyplot
        """
        # TODO: add functionality for drawing a horizontal line where a Qini
        # curve plateaus
        try:
            import matplotlib.pyplot as plt
        except:
            raise ImportError("plot method requires matplotlib.")

        if not "color" in kwargs:
            kwargs["color"] = "black"
        plt.plot(self.path_spend_, self.path_gain_, **kwargs)
        if "label" in kwargs:
            plt.legend(loc="upper left")
        if show_ci:
            ub = self.path_gain_ + 1.96 * self.path_std_err_
            lb = self.path_gain_ - 1.96 * self.path_std_err_
            plt.plot(
                self.path_spend_,
                ub,
                color=kwargs["color"],
                linestyle="dotted",
                linewidth=1,
            )
            plt.plot(
                self.path_spend_,
                lb,
                color=kwargs["color"],
                linestyle="dotted",
                linewidth=1,
            )
        plt.xlabel("spend")
        plt.ylabel("gain")

    @property
    def path_spend_(self):
        self._ensure_fit()
        return self._path["spend"]

    @property
    def path_gain_(self):
        self._ensure_fit()
        return self._path["gain"]

    @property
    def path_std_err_(self):
        self._ensure_fit()
        return self._path["std_err"]

    @property
    def path_allocated_unit_(self):
        self._ensure_fit()
        return self._path["ipath"]

    @property
    def path_allocated_arm_(self):
        self._ensure_fit()
        return self._path["kpath"]

    def _ensure_fit(self):
        if not self._is_fit:
            raise ValueError("MAQ object is not fit.")

    def __repr__(self):
        if self._is_fit:
            return "MAQ object with {} units and {} arms.".format(
                self._dim[0], self._dim[1]
            )
        else:
            return "MAQ object (not fit)."
