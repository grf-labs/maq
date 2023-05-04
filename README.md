# maq

[![Build Status](https://dev.azure.com/grf-labs/grf/_apis/build/status/grf-labs.maq?branchName=master)](https://dev.azure.com/grf-labs/grf/_build/latest?definitionId=5&branchName=master)

A package for evaluating treatment targeting with multiple arms via the Multi-Action Qini ("maq") - a generalization of the Qini to multiple costly treatment arms.

The development version (R package) can be installed by

```R
devtools::install_github("grf-labs/maq", subdir = "r-package/maq")
```
(Installing from source requires a compiler that implements C++11 or later).

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
cost.hat <- X[eval, 4:5]

# Fit an evaluation forest to compute doubly robust evaluation set scores.
eval.forest <- grf::multi_arm_causal_forest(X[eval, ], Y[eval], W[eval])
DR.scores <- grf::get_scores(eval.forest)[,,]

# Fit a MAQ using evaluation set estimates.
max.budget <- 1
mq <- maq(tau.hat, cost.hat, max.budget, DR.scores)

# Plot the MAQ curve.
plot(mq)

# Get an estimate of optimal reward at a given spend per unit along with standard errors.
average_gain(mq, spend = 0.3)

# Get the optimal treatment allocation matrix at a given spend per unit.
pi.mat <- predict(mq, spend = 0.3)

# If the treatment randomization probabilities are known, then an alternative to
# evaluation via AIPW scores is to use inverse-propensity weighting (IPW).
W.hat.true <- rep(1/3, 3)
observed.W <- match(W, levels(W))
Y.k.mat <- matrix(0, length(W), nlevels(W))
Y.k.mat[cbind(seq_along(observed.W), observed.W)] <- Y
Y.k.ipw <- sweep(Y.k.mat, 2, W.hat.true, "/")
Y.k.ipw.eval <- Y.k.ipw[eval, -1] - Y.k.ipw[eval, 1]

mq.ipw <- maq(tau.hat, cost.hat, max.budget, Y.k.ipw.eval)
```

### Details

Consider a set of costly and mutually exclusive treatment arms $k = 0, \ldots, K$ where $k=0$ is a zero-cost control. Let $\tau(X_i)$ be a vector of treatment effects for unit $i$, i.e. the $k$-th element ($k > 0$) is $\mu_{ik} - \mu_{i0}$, where $\mu_{ik} = E[Y_i(k) | X_i = x]$. Let $C(X_i)$ be a vector of positive costs, i.e. the $k$-th element is the cost of assigning unit $i$ arm $k$.

The multi-action Qini is defined as the value of the optimal cost-constrained treatment allocation $\pi_b(X_i) \in [0, 1]^K$ at any budget constraint $b$. `maq` delivers the path of these optimal allocations and confidence intervals for the optimal value (on a held out evaluation set) by efficiently solving the following series of linear programs:

```math
\begin{aligned}
\max_{\pi_b} \quad & \mathbb{E}[\langle \pi_b(X_i),~ \tau(X_i) \rangle] \\
\textrm{s.t.} \quad & \mathbb{E}[\langle \pi_b(X_i),~ C(X_i) \rangle] \leq b \\
& \langle \pi_b(X_i),~ \mathbf{1} \rangle \leq 1 \\
& \pi_b(X_i) \geq 0.
\end{aligned}
```

### References
