library(microbenchmark)
library(maq)
rm(list = ls())
set.seed(42)

# 1 medium n medium K
n = 1e5
K = 10
reward = matrix(1 + rnorm(n*K), n, K)
cost = 0.05 + matrix(runif(n*K), n, K)

b1 <- microbenchmark(
  maq(reward, cost, reward, R = 200, paired.inference = FALSE),
  times = 10,
  unit = "seconds"
)

# 2 medium n large K
n = 1e5
K = 100
reward = matrix(1 + rnorm(n*K), n, K)
cost = 0.05 + matrix(runif(n*K), n, K)

b2 <- microbenchmark(
  maq(reward, cost, reward, R = 200, paired.inference = FALSE),
  times = 10,
  unit = "seconds"
)

# 3 medium n extreme K
n = 1e5
K = 1000
reward = matrix(1 + rnorm(n*K), n, K)
cost = 0.05 + matrix(runif(n*K), n, K)

b3 <- microbenchmark(
  maq(reward, cost, reward, R = 200, paired.inference = FALSE),
  times = 5,
  unit = "seconds"
)

# 4 large n medium K
n = 1e6
K = 10
reward = matrix(1 + rnorm(n*K), n, K)
cost = 0.05 + matrix(runif(n*K), n, K)

b41 <- microbenchmark(
  maq(reward, cost, reward, R = 0, paired.inference = FALSE),
  times = 10,
  unit = "seconds"
)

b42 <- microbenchmark(
  maq(reward, cost, reward, R = 200, paired.inference = FALSE),
  times = 10,
  unit = "seconds"
)

# 5 large n large K
n = 1e6
K = 100
reward = matrix(1 + rnorm(n*K), n, K)
cost = 0.05 + matrix(runif(n*K), n, K)

b5 <- microbenchmark(
  maq(reward, cost, reward, R = 0, paired.inference = FALSE),
  times = 5,
  unit = "seconds"
)

print(b1, digits = 3)
print(b2, digits = 3)
print(b3, digits = 3)
print(b41, digits = 3)
print(b42, digits = 3)
print(b5, digits = 3)
print(list(b1,b2,b3,b41,b42,b5), digits = 3)
