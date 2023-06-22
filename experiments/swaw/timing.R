rm(list = ls())
set.seed(42)
library(maq)
library(microbenchmark)

n = 1e6
K = 5
reward = matrix(runif(n * K), n, K)
reward.eval = matrix(runif(n * K), n, K)
cost = 0.05 + matrix(runif(n * K), n, K)

print(microbenchmark(
  mq <- maq(reward, cost, 1e5, cost, R = 0),
  times = 100,
  unit = "seconds"
), digits = 2)
# Unit: seconds
#                                         expr min  lq mean median  uq max neval
# mq <- maq(reward, cost, 1e+05, cost, R = 0) 1.4 1.4  1.5    1.5 1.6 1.6   100
