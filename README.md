# Multi-Armed Qini

[![CRANstatus](https://www.r-pkg.org/badges/version/maq)](https://cran.r-project.org/package=maq)
[![Build Status](https://dev.azure.com/grf-labs/grf/_apis/build/status/grf-labs.maq?branchName=master)](https://dev.azure.com/grf-labs/grf/_build/latest?definitionId=5&branchName=master)

A package for policy evaluation using generalized Qini curves: Evaluate data-driven treatment targeting rules for one or more treatment arms over different budget constraints in experimental or observational settings under unconfoundedness.

* Introduction: [Qini curves: Automatic cost-benefit analysis.](https://grf-labs.github.io/grf/articles/maq.html)

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
ma.qini <- maq(tau.hat, cost, DR.scores, R = 200)

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
mq.ipw <- maq(tau.hat, cost, IPW.scores)
```

### Details

Let $\hat \tau(\cdot)$ be an estimated CATE function where the $k$-th element measures the conditional average treatment effect $E[Y_i(k) - Y_i(0) | X_i]$ for a given unit $X_i$ for one of $k=1, \ldots, K$ treatment arms, where $k=0$ is a control arm. Let $C(\cdot)$ be some known cost function that quantifies the cost of assigning a given treatment arm to the $i$-th unit. `maq` delivers estimates of the Qini curve

$$
Q(B) = E[\langle \pi_B(X_i),~ \tau(X_i)\rangle],
$$

which is the expected gain, at any budget constraint $B$, when assigning treatment using the policy $\pi_B$ that optimally selects (using the given functions $\hat \tau(\cdot)$ and $C(\cdot))$ which arm to assign to which unit such that the average incurred cost is less than or equal to $B$. The policy $\pi_B$ is a solution to a linear program: `maq` computes a solution path for these treatment allocations over increasing budget levels $B$ via an algorithm that leverages the multiple-choice knapsack structure of this problem. See [the algorithm reference](https://github.com/grf-labs/maq/tree/master/REFERENCE.md) for an overview.

### References

Erik Sverdrup, Han Wu, Susan Athey, and Stefan Wager.
<b>Qini Curves for Multi-Armed Treatment Rules.</b>
<i>Journal of Computational and Graphical Statistics</i>, forthcoming.
[<a href="https://doi.org/10.1080/10618600.2024.2418820">paper</a>,
<a href="https://arxiv.org/abs/2306.11979">arxiv</a>]
