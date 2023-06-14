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
budget = 1e9
n = 1e5
K = 10
reward = matrix(1 + rnorm(n*K), n, K)
cost = 0.05 + matrix(runif(n*K), n, K)

b1 <- microbenchmark(
  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE),
  times = 10,
  unit = "seconds"
)

# 2 medium n large K
n = 1e5
K = 100
reward = matrix(1 + rnorm(n*K), n, K)
cost = 0.05 + matrix(runif(n*K), n, K)

b2 <- microbenchmark(
  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE),
  times = 10,
  unit = "seconds"
)

# 3 medium n extreme K
n = 1e5
K = 1000
reward = matrix(1 + rnorm(n*K), n, K)
cost = 0.05 + matrix(runif(n*K), n, K)

b3 <- microbenchmark(
  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE),
  times = 5,
  unit = "seconds"
)

# 4 large n medium K
n = 1e6
K = 10
reward = matrix(1 + rnorm(n*K), n, K)
cost = 0.05 + matrix(runif(n*K), n, K)

b41 <- microbenchmark(
  maq(reward, cost, budget, reward, R = 0, paired.inference = FALSE),
  times = 10,
  unit = "seconds"
)

b42 <- microbenchmark(
  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE),
  times = 10,
  unit = "seconds"
)

# 5 large n large K
n = 1e6
K = 100
reward = matrix(1 + rnorm(n*K), n, K)
cost = 0.05 + matrix(runif(n*K), n, K)

b5 <- microbenchmark(
  maq(reward, cost, budget, reward, R = 0, paired.inference = FALSE),
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


#' @ Do weights/tie at runtime (#15) https://github.com/grf-labs/maq/pull/15
# [[1]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 2.61 2.71 3.01   3.04 3.23 3.38    10
#
# [[2]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 4.42 4.72  5.3   5.27 5.74 6.22    10
#
# [[3]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 17.4 17.6 17.8   17.7 17.9 18.2     5
#
# [[4]]
# Unit: seconds
# expr  min   lq mean median   uq max neval
# maq(reward, cost, budget, reward, R = 0, paired.inference = FALSE) 1.96 1.96 2.18   2.02 2.42 2.8    10
#
# [[5]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 31.5 31.5 32.1     32 32.3 33.2    10
#
# [[6]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 0, paired.inference = FALSE) 10.8 10.8 10.9   10.8 10.9 11.1     5


#' @ Use plain std::sort
# [[1]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 1.97 2.02 2.03   2.04 2.05 2.08    10
#
# [[2]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 3.79 3.8 3.81   3.81 3.82 3.83    10
#
# [[3]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 15.2 15.3 15.3   15.3 15.4 15.4     5
#
# [[4]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 0, paired.inference = FALSE) 1.96 1.97 1.99   1.98 1.99 2.05    10
#
# [[5]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 31.4 31.6 31.8   31.7 32.1 32.4    10
#
# [[6]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 0, paired.inference = FALSE) 10.2 10.2 10.3   10.2 10.5 10.5     5


#' @ Add templates (#10) https://github.com/grf-labs/maq/pull/12
# [[1]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 2.09 2.11 2.17   2.16 2.2 2.3    10
#
# [[2]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 4.82 4.85 4.89   4.88 4.92   5    10
#
# [[3]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 19.9 19.9   20   19.9 20 20.1     5
#
# [[4]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 0, paired.inference = FALSE) 1.93 1.94 1.97   1.96 1.98 2.05    10
#
# [[5]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 31.8 31.9 32.1     32 32.3 32.7    10
#
# [[6]]
# Unit: seconds
#  maq(reward, cost, budget, reward, R = 0, paired.inference = FALSE) 16.6 16.6 16.7   16.7 16.7 16.7     5


#' @ add pibar (#10) https://github.com/grf-labs/maq/pull/10
# [[1]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 1.98 2.03 2.06   2.06 2.07 2.13    10
#
# [[2]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 5.86 5.89 5.94   5.94 5.95 6.08    10
#
# [[3]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 35.2 35.3 35.3   35.3 35.4 35.4     5
#
# [[4]]
# Unit: seconds
# expr min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 0, paired.inference = FALSE) 2.3 2.31 2.34   2.33 2.34 2.42    10
#
# [[5]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200, paired.inference = FALSE) 31.8 31.9 32.4   32.3 32.6 33.4    10
#
# [[6]]
# Unit: seconds
# expr min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 0, paired.inference = FALSE)  31 31.1   38   31.1 32.1 64.8     5


#' @ 9edff55064ca75b7806b0532ec24f82096b22010
# [[1]]
# Unit: seconds
# expr min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200) 1.9 1.91 1.93   1.92 1.94 1.96    10
#
# [[2]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200) 4.64 4.65 4.67   4.67 4.69 4.72    10
#
# [[3]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200) 22.2 22.2 22.3   22.3 22.3 22.5     5
#
# [[4]]
# Unit: seconds
# expr min lq mean median   uq max neval
# maq(reward, cost, budget, reward, R = 0)   2  2 2.02   2.01 2.02 2.1    10
#
# [[5]]
# Unit: seconds
# expr  min lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 200) 31.7 32 32.1     32 32.1 32.4    10
#
# [[6]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, reward, R = 0) 18.2 18.2 18.2   18.3 18.3 18.3     5


#' @ ec1d19cb19b9001886ba468e638edaece4b40579
# [[1]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, R = 200) 1.81 1.83 1.85   1.84 1.89 1.91    10
#
# [[2]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, R = 200) 4.32 4.32 4.34   4.33 4.35 4.36    10
#
# [[3]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, R = 200) 21.4 21.5 21.5   21.5 21.6 21.7     5
#
# [[4]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, R = 0) 1.95 1.98 1.99   1.99 2.01 2.04    10
#
# [[5]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, R = 200) 29.5 29.6 29.9   29.7 30.1 31.4    10
#
# [[6]]
# Unit: seconds
# expr  min   lq mean median   uq  max neval
# maq(reward, cost, budget, R = 0) 18.1 18.1 18.2   18.1 18.2 18.3     5
