# maq
Python bindings for `maq` can be installed from source by

```
pip install git+https://github.com/grf-labs/maq.git#egg=maq\&subdirectory=python-package
```

(requires Cython/Numpy)

### Example

```python
import numpy as np
from maq import MAQ

n = 500
K = 4
reward = np.random.randn(n, K)
cost = np.random.rand(n, K)

# Fit a MAQ up to the maximum spend per unit.
mq = MAQ(n_bootstrap=150)
mq.fit(reward, cost, np.mean(cost))

# Get an estimate of optimal reward along with standard errors.
mq.average_gain(spend=0.1)

# Get the optimal treatment allocation matrix at a given spend.
mq.predict(spend=0.1)
```
