# Check that the following code paths passes under valgrind.
# Usage:
# R -d "valgrind --tool=memcheck --leak-check=full" --vanilla  < test_maq_valgrind.R
library(maq)
budget <- 1000
n <- 1500
K <- 10
reward <- matrix(0.1 + rnorm(n * K), n, K)
reward.eval <- matrix(0.1 + rnorm(n * K), n, K)
cost <- 0.05 + matrix(runif(n * K), n, K)
wts <- runif(n)
clust <- sample(1:250, n, TRUE)

mq <- maq(reward, cost, budget, reward.eval, sample.weights = wts, clusters = clust, R = 200)

mqr <- maq(reward, cost, budget, reward.eval, sample.weights = wts, clusters = clust, R = 200, target.with.covariates = FALSE)
