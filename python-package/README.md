# maq

### Installation

Python bindings for `maq` can be installed from source with

```
pip install "git+https://github.com/grf-labs/maq.git#egg=maq&subdirectory=python-package"
```

Compiling from source requires Cython/Numpy. There are no C++ dependencies except the Standard Template Library, and C++11 or higher.

Currently, only the core MAQ solver features (point estimates and confidence intervals of the solution path) are implemented. The R package has support for all the functionality (like clustered/paired bootstraps, etc).

### Usage Example

```python
import numpy as np
from maq import MAQ

n = 500
K = 4
reward = np.random.randn(n, K)
cost = np.random.rand(n, K)
reward_eval = np.random.randn(n, K)

# Fit a MAQ up to the maximum spend per unit.
max_budget = np.mean(cost)
mq = MAQ(budget=max_budget, n_bootstrap=200)
mq.fit(reward, cost, reward_eval)

# Get an estimate of optimal reward along with standard errors.
mq.average_gain(spend=0.1)

# Get the optimal treatment allocation matrix at a given spend.
mq.predict(spend=0.1)

# Plot the gain curve.
import matplotlib.pyplot as plt

plt.plot(mq.path_spend_, mq.path_gain_)
plt.xlabel("Spend/unit")
plt.title("Gain/unit")
plt.show()

# Show 95% confidence bars.
ub = mq.path_gain_ + 1.96 * mq.path_std_err_
lb = mq.path_gain_ - 1.96 * mq.path_std_err_
plt.plot(mq.path_spend_, ub, color="black", linestyle="dashed")
plt.plot(mq.path_spend_, lb, color="black", linestyle="dashed")
plt.show()
```
