# maq

[![Build Status](https://dev.azure.com/grf-labs/grf/_apis/build/status/grf-labs.maq?branchName=master)](https://dev.azure.com/grf-labs/grf/_build/latest?definitionId=5&branchName=master)

A C++ solver for the Multi-Action QINI (maq).

Consider a set of costly and mutually exclusive treatment arms $k = 0, \ldots, K - 1$ where $k=0$ is a zero-cost control. Let $\tau(X_i)$ be a vector of treatment effects for unit $i$, i.e. the $k$-th element is $\mu_{ik} - \mu_{i0}$, where $\mu_{ik} = E[Y_i(k) | X_i = x]$. Let $C(X_i)$ be a vector of costs, i.e. the $k$-th element is the cost of assigning user $i$ arm $k$. `maq` finds the optimal cost-constrained treatment allocation $\pi_b(X_i) \in [0, 1]^{K-1}$ for any budget constraint $b$.

In particular, `maq` efficiently computes the path of solutions to the following stochastic LP:

```math
\begin{aligned}
\max_{\pi_b} \quad & \mathbb{E}[\langle \pi_b(X_i),~ \tau(X_i) \rangle] \\
\textrm{s.t.} \quad & \mathbb{E}[\langle \pi_b(X_i),~ C(X_i) \rangle] \leq b \\
& \langle \pi_b(X_i),~ \mathbf{1} \rangle \leq 1 \\
& \pi_b(X_i) \geq 0,
\end{aligned}
```

by using a priority-queue based algorithm that leverages the Knapsack-type structure of this problem.

The development version (R) can be installed by

```r
devtools::install_github("grf-labs/maq", subdir = "r-package/maq")
```

Python bindings are [here](https://github.com/grf-labs/maq/tree/master/python-package).

### Usage Example

### References
