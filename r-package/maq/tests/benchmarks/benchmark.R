library(microbenchmark)
library(maq)
rm(list = ls())
set.seed(42)

sessionInfo()
# R version 4.2.2 (2022-10-31)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.6.3
parallel::detectCores()
# [1] 8

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

# v0.1.0
# > print(list(b1,b2,b3,b41,b42,b5), digits = 3)
# [[1]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 200, paired.inference = FALSE) 1.93 1.97    2   1.99 2.02 2.12    10
#
# [[2]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 200, paired.inference = FALSE) 3.84 3.85 3.87   3.87 3.89 3.92    10
#
# [[3]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 200, paired.inference = FALSE) 15.3 15.3 15.4   15.3 15.3 15.8     5
#
# [[4]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 0, paired.inference = FALSE) 1.91 1.92 1.94   1.93 1.95 2.01    10
#
# [[5]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 200, paired.inference = FALSE) 31.7 31.8 32.9   32.1 32.4 40.7    10
#
# [[6]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, reward, R = 0, paired.inference = FALSE) 10.2 10.2 10.2   10.2 10.2 10.2     5
