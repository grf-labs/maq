# Solves the linear knapsack-type problem with a generic LP solver for a given spend.
# attachNamespace("lpSolve")
lp_solver = function(reward, cost, budget) {
  library(lpSolve)
  K = ncol(reward)
  x.coeffs = c(t(reward)) / NROW(reward)
  A.mat = matrix(0, nrow(reward), length(x.coeffs))
  start = 1
  for (row in 1:nrow(A.mat)) {
    end = start + K -1
    A.mat[row, start:end] = 1
    start = end + 1
  }
  c.coeff = c(t(cost)) / NROW(cost)
  f.con = rbind(A.mat, c.coeff)
  f.dir = rep("<=", nrow(f.con))
  f.rhs = c(rep(1, nrow(A.mat)), budget)
  lp.res = lp("max", x.coeffs, f.con, f.dir, f.rhs)

  list(gain = sum(x.coeffs * lp.res$solution),
       spend = sum(c.coeff * lp.res$solution),
       alloc.mat = matrix(lp.res$solution, nrow = nrow(reward), byrow = TRUE))
}

test_that("gain solution path is same as LP solution (small problem)", {
  budget <- 1e6
  n <- 50
  K <- 10
  reward <- matrix(1 + 2 * runif(n * K), n, K);
  cost <- matrix(runif(n * K), n, K);

  mqini <- maq(reward, cost, budget)
  gain <- mqini[["_path"]]$gain

  lp.gain <- c()
  for (ix in seq_along(mqini[["_path"]]$spend)) {
    spend <- mqini[["_path"]]$spend[ix]
    lp.gain <- c(lp.gain, lp_solver(reward, cost, spend)$gain)
  }

  expect_equal(gain, lp.gain, tolerance = 1e-12)
})

test_that("solution is same as LP solution (fixed medium problem)", {
  budget <- 5
  n <- 500
  K <- 10
  reward <- matrix(1 + 2 * runif(n * K), n, K);
  cost <- matrix(runif(n * K), n, K);

  mqini <- maq(reward, cost, budget)

  ix <- length(mqini[["_path"]]$spend) - 1
  spend <- mqini[["_path"]]$spend[ix]
  gain <- mqini[["_path"]]$gain[ix]

  expect_equal(gain, !!lp_solver(reward, cost, spend)$gain, tolerance = 1e-12)
})

test_that("pi matrix is consistent with gain path", {
  budget <- 25
  n <- 100
  K <- 10
  reward <- matrix(1 + 2 * runif(n * K), n, K);
  cost <- matrix(runif(n * K), n, K);

  mqini <- maq(reward, cost, budget)
  gain <- mqini[["_path"]]$gain

  pi.gain <- c()
  for (ix in seq_along(mqini[["_path"]]$spend)) {
    spend <- mqini[["_path"]]$spend[ix]
    pi.mat <- predict(mqini, spend = spend)
    pi.gain <- c(pi.gain, sum(reward * pi.mat) / n)
  }
  expect_equal(pi.gain, gain, tolerance = 1e-12)
})

test_that("arbitrary points off gain path is same as LP", {
  budget <- 1e6
  n <- 500
  K <- 5
  reward <- matrix(0.1 + rnorm(n * K), n, K);
  cost <- matrix(0.05 + runif(n * K), n, K);

  mqini <- maq(reward, cost, budget)
  spend.rng <- range(mqini[["_path"]]$spend)
  spends <- runif(25, spend.rng[1], spend.rng[2])

  gain.mq <- c()
  gain.mq.alt <- c()
  gain.lp <- c()
  for (spend in spends) {
    gain.est <- average_gain(mqini, spend = spend)[[1]]
    pi.mat <- predict(mqini, spend = spend)
    gain.est.alt <- sum(reward * pi.mat) / n
    gain.est.lp <- lp_solver(reward, cost, spend)$gain

    gain.mq <- c(gain.mq, gain.est)
    gain.mq.alt <- c(gain.mq.alt, gain.est.alt)
    gain.lp <- c(gain.lp, gain.est.lp)
  }

  expect_equal(gain.mq, gain.lp, tolerance = 1e-12)
  expect_equal(gain.mq.alt, gain.lp, tolerance = 1e-12)
})
