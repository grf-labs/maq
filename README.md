# maq

[![Build Status](https://dev.azure.com/grf-labs/grf/_apis/build/status/grf-labs.maq?branchName=master)](https://dev.azure.com/grf-labs/grf/_build/latest?definitionId=5&branchName=master)

A C++ solver for the Multi-Action QINI ("maq") - a generalization of the QINI to multiple costly treatment arms.

The development version (R package) can be installed by

```R
devtools::install_github("grf-labs/maq", subdir = "r-package/maq")
```

Python bindings are [here](https://github.com/grf-labs/maq/tree/master/python-package).

### Usage Example

```R
# Fit a CATE estimator (using GRF) on a training sample.
n <- 2000
p <- 5
X <- matrix(runif(n * p), n, p)
W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
Y <- X[, 1] + X[, 2] * (W == "B") + X[, 3] * (W == "C") + rnorm(n)
train <- sample(1:n, n/2)
eval <- -train

tau.forest <- grf::multi_arm_causal_forest(X[train, ], Y[train], W[train])

# Predict CATEs on held out evaluation data.
tau.hat <- predict(tau.forest, X[eval, ])$predictions[,,]

# Form cost estimates - the following are a toy example.
cost.hat <- X[eval, 2:3]

# Fit a evaluation forest to compute doubly robust evaluation set scores.
eval.forest <- grf::multi_arm_causal_forest(X[eval, ], Y[eval], W[eval])
DR.scores <- grf::get_scores(eval.forest)[,,]

# Fit a MAQ using evaluation set estimates.
max.budget <- 1
mq <- maq(tau.hat, cost.hat, max.budget, DR.scores)

# Plot the MAQ curve.
plot(mq)

# Get an estimate of optimal reward at a given spend/unit along with standard errors.
average_gain(mq, spend = 0.3)

# Get the optimal treatment allocation matrix at a given spend/unit.
pi.mat <- predict(mq, spend = 0.3)
```

### Details

Consider a set of costly and mutually exclusive treatment arms $k = 0, \ldots, K$ where $k=0$ is a zero-cost control. Let $\tau(X_i)$ be a vector of treatment effects for unit $i$, i.e. the $k$-th element ($k > 0$) is $\mu_{ik} - \mu_{i0}$, where $\mu_{ik} = E[Y_i(k) | X_i = x]$. Let $C(X_i)$ be a vector of positive costs, i.e. the $k$-th element is the cost of assigning user $i$ arm $k$. `maq` finds the optimal cost-constrained treatment allocation $\pi_b(X_i) \in [0, 1]^K$ for any budget constraint $b$.

In particular, `maq` efficiently computes the path of solutions to the following stochastic LP:

```math
\begin{aligned}
\max_{\pi_b} \quad & \mathbb{E}[\langle \pi_b(X_i),~ \tau(X_i) \rangle] \\
\textrm{s.t.} \quad & \mathbb{E}[\langle \pi_b(X_i),~ C(X_i) \rangle] \leq b \\
& \langle \pi_b(X_i),~ \mathbf{1} \rangle \leq 1 \\
& \pi_b(X_i) \geq 0.
\end{aligned}
```

### References
