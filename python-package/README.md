# maq
Python bindings for `maq` can be installed from source by

```
pip install "git+https://github.com/grf-labs/maq.git#egg=maq&subdirectory=python-package"
```

(requires Cython/Numpy)

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
mq = MAQ(n_bootstrap=200)
max_budget = np.mean(cost)
mq.fit(reward, cost, max_budget, reward_eval)

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
