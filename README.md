# maq

[![CRANstatus](https://www.r-pkg.org/badges/version/maq)](https://cran.r-project.org/package=maq)
[![Build Status](https://dev.azure.com/grf-labs/grf/_apis/build/status/grf-labs.maq?branchName=master)](https://dev.azure.com/grf-labs/grf/_build/latest?definitionId=5&branchName=master)

A package for evaluating multi-armed treatment rules via the Multi-Armed Qini ("maq") - a generalization of the Qini to multiple costly treatment arms.

* Introduction: [Qini curves: Automatic cost-benefit analysis.](https://grf-labs.github.io/grf/articles/maq.html)

* [Changelog](https://github.com/grf-labs/maq/tree/master/CHANGELOG.md).

### Installation

The latest release of the package can be installed through CRAN:

```R
install.packages("maq")
```

The development version can be installed with:

```R
devtools::install_github("grf-labs/maq", subdir = "r-package/maq")
```
(Installing from source requires a compiler that implements C++11 or later)

**Python** bindings are [here](https://github.com/grf-labs/maq/tree/master/python-package).

### Usage Example

```R
library(maq)

# Fit a CATE estimator on a training sample.
n <- 3000
p <- 5
X <- matrix(runif(n * p), n, p)
W <- as.factor(sample(c("0", "1", "2"), n, replace = TRUE))
Y <- X[, 1] + X[, 2] * (W == "1") + 1.5 * X[, 3] * (W == "2") + rnorm(n)
train <- sample(1:n, n/2)

tau.forest <- grf::multi_arm_causal_forest(X[train, ], Y[train], W[train])

# Predict CATEs on held out evaluation data.
test <- -train
tau.hat <- predict(tau.forest, X[test, ], drop = TRUE)$predictions

# Assume costs equal a unit's pre-treatment covariate - the following are a toy example.
cost <- cbind(X[test, 4] / 4, X[test, 5])

# Fit an evaluation forest to compute doubly robust scores on the test set.
eval.forest <- grf::multi_arm_causal_forest(X[test, ], Y[test], W[test])
DR.scores <- grf::get_scores(eval.forest, drop = TRUE)

# Fit a Qini curve on evaluation data, using 200 bootstrap replicates for confidence intervals.
max.budget <- 1
ma.qini <- maq(tau.hat, cost, max.budget, DR.scores, R = 200)

# Plot the Qini curve.
plot(ma.qini)
legend("topleft", c("All arms", "95% CI"), lty = c(1, 3))

# Get an estimate of gain at a given spend per unit along with standard errors.
average_gain(ma.qini, spend = 0.2)

# Get the treatment allocation matrix at a given spend per unit.
pi.mat <- predict(ma.qini, spend = 0.2)

# If the treatment randomization probabilities are known, then an alternative to
# evaluation via AIPW scores is to use inverse-propensity weighting (IPW).
W.hat <- rep(1/3, 3)
observed.W <- match(W, levels(W))
Y.mat <- matrix(0, length(W), nlevels(W))
Y.mat[cbind(seq_along(observed.W), observed.W)] <- Y
Y.ipw <- sweep(Y.mat, 2, W.hat, "/")
Y.ipw.test <- Y.ipw[test, -1] - Y.ipw[test, 1]

mq.ipw <- maq(tau.hat, cost, max.budget, Y.ipw.test)
```

### Details

Consider a set of costly and mutually exclusive treatment arms $k = 0, \ldots, K$ where $k=0$ is a zero-cost control. Let $\hat \tau(X_i)$ be a vector of treatment effects estimates for unit $i$, i.e. the $k$-th element ($k > 0$) is $\hat \mu_{ik} - \hat \mu_{i0}$, where $\mu_{ik} = E[Y_i(k) | X_i = x]$. Let $\widehat C(X_i)$ be a vector of positive cost estimates, i.e. the $k$-th element is the estimated cost of assigning unit $i$ arm $k$.

The multi-armed Qini is defined as the value of the optimal cost-constrained treatment allocation $\pi_B(X_i) \in [0, 1]^K$ at any budget constraint $B$. `maq` delivers the path of these optimal allocations and confidence intervals for the estimated value (on a held-out test set) by efficiently solving a series of linear programs, for each $B \in (0, B_{max}]$:

```math
\begin{aligned}
\max_{\pi_B} \quad & \sum_{i=1}^{n} \langle \pi_B(X_i),~ \hat \tau(X_i) \rangle \\
\textrm{s.t.} \quad & \sum_{i=1}^{n} \langle \pi_B(X_i),~ \widehat C(X_i) \rangle \leq B \\
& \langle \pi_B(X_i),~ \mathbf{1} \rangle \leq 1 \\
& \pi_B(X_i) \geq 0.
\end{aligned}
```

### References

Erik Sverdrup, Han Wu, Susan Athey, and Stefan Wager.
<b>Qini Curves for Multi-Armed Treatment Rules.</b> 2023.
[<a href="https://arxiv.org/abs/2306.11979">arxiv</a>]
