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
IPW.scores <- get_ipw_scores(Y[test], W[test], W.hat)
mq.ipw <- maq(tau.hat, cost, max.budget, IPW.scores)
```

### Details

Consider a set of costly and mutually exclusive treatment arms $k = 0, \ldots, K$ where $k=0$ is a zero-cost control. Let $\hat \tau(X_i)$ be a vector of treatment effect estimates for unit $i$, i.e. the $k$-th element ($k > 0$) is $\hat \mu_{ik} - \hat \mu_{i0}$, where $\mu_{ik} = E[Y_i(k) | X_i = x]$. Let $C(X_i)$ be a vector of positive cost, i.e. the $k$-th element measures the cost of assigning unit $i$ arm $k$.

Given the functions $\hat \tau(X_i)$ and $C(\cdot)$, `maq` delivers a Qini curve that quantifies the value of optimally assigning treatment while satisfying a budget constraint $B$. The policies underlying this treatment allocation are a solution to a series of linear programs, which `maq` solves efficiently.

### References

Erik Sverdrup, Han Wu, Susan Athey, and Stefan Wager.
<b>Qini Curves for Multi-Armed Treatment Rules.</b> 2023.
[<a href="https://arxiv.org/abs/2306.11979">arxiv</a>]
