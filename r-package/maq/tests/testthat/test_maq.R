test_that("maq works as expected", {
  budget <- 1000
  n <- 500
  K <- 3
  reward <- matrix(0.1 + rnorm(n * K), n, K)
  cost <- 0.05 + matrix(runif(n * K), n, K)

  mq <- maq(-abs(reward), cost, budget, R = 150)
  plot(mq)
  expect_true(all(predict(mq, 10) == 0))
  expect_equal(average_gain(mq, 10), c(estimate = 0, std.err = 0))

  mq <- maq(reward[, 1], cost[, 1], budget, R = 150)
  Matrix::which(predict(mq, 10) > 0)

  # scale invariances
  mq <- maq(reward, cost, budget, seed = 42)
  mq.scale <- maq(reward * 1000, cost * 1000, budget, seed = 42)

  expect_equal(mq[["_path"]]$spend, mq.scale[["_path"]]$spend / 1000, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$gain, mq.scale[["_path"]]$gain / 1000, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$std.err, mq.scale[["_path"]]$std.err / 1000, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$ipath, mq.scale[["_path"]]$ipath)
  expect_equal(mq[["_path"]]$kpath, mq.scale[["_path"]]$kpath)
})

test_that("sample weighting works as expected", {
  budget <- 1000
  n <- 200
  K <- 3
  reward <- matrix(0.1 + rnorm(n * K), n, K)
  cost <- 0.05 + matrix(runif(n * K), n, K)

  # giving weight 2 ~same sample as duplicating
  dupe <- sample(1:n, 100)
  reward.dupe <- rbind(reward, reward[dupe, ])
  cost.dupe <- rbind(cost, cost[dupe, ])
  wts <- rep(1, n)
  wts[dupe] <- 2

  mq <- maq(reward, cost, budget, sample.weights = wts)
  mq.dupe <- maq(reward.dupe, cost.dupe, budget)

  spends <- c(0.1, 0.25, 0.3, 0.35, 0.4, 0.5)
  est <- lapply(spends, function(s) average_gain(mq, s)[[1]])
  est.dupe <- lapply(spends, function(s) average_gain(mq.dupe, s)[[1]])

  expect_equal(est, est.dupe, tolerance = 0.05)
})

test_that("clustering works as expected", {
  budget <- 1000
  n <- 100
  K <- 3
  reward <- matrix(0.1 + rnorm(n * K), n, K)
  cost <- 0.05 + matrix(runif(n * K), n, K)

  rewardc <- rbind(reward, reward, reward, reward, reward)
  costc <- rbind(cost, cost, cost, cost, cost)
  clust <- rep(1:n, 5)

  mq <- maq(reward, cost, budget)
  mq.clust <- maq(rewardc, costc, budget, clusters = clust)

  spends <- c(0.1, 0.25, 0.3, 0.35, 0.4, 0.5)
  est <- lapply(spends, function(s) average_gain(mq, s))
  est.clust <- lapply(spends, function(s) average_gain(mq.clust, s))

  expect_equal(est, est.clust, tolerance = 0.075)
})

test_that("std errors works as expected", {
  # Get exact coverage of points on the curve
  budget <- 100
  ntrue <- 100000
  K <- 5
  reward <- matrix(0.1 + rnorm(ntrue * K), ntrue, K)
  cost <- 0.05 + matrix(runif(ntrue * K), ntrue, K)
  mqt <- maq(reward, cost, budget, R = 0)

  spend1 <- 0.1
  spend2 <- 0.25
  spend3 <- 0.4
  true1 <- average_gain(mqt, spend = spend1)[[1]]
  true2 <- average_gain(mqt, spend = spend2)[[1]]
  true3 <- average_gain(mqt, spend = spend3)[[1]]

  n <- 1000
  res <- t(replicate(500, {
    reward <- matrix(0.1 + rnorm(n * K), n, K)
    cost <- 0.05 + matrix(runif(n * K), n, K)
    mq <- maq(reward, cost, budget, R = 150)

    pp1 <- average_gain(mq, spend1)
    pp2 <- average_gain(mq, spend2)
    pp3 <- average_gain(mq, spend3)

    c(
      est1 = pp1[[1]],
      se1 = pp1[[2]],
      bias1 = pp1[[1]] - true1,
      cov1 = abs(pp1[[1]] - true1) / pp1[[2]] <= 1.96,

      est2 = pp2[[1]],
      se2 = pp2[[2]],
      bias2 = pp2[[1]] - true2,
      cov2 = abs(pp2[[1]] - true2) / pp2[[2]] <= 1.96,

      est3 = pp3[[1]],
      se3 = pp3[[2]],
      bias3 = pp3[[1]] - true3,
      cov3 = abs(pp3[[1]] - true3) / pp3[[2]] <= 1.96
    )
  }))

  expect_equal(mean(res[, "cov1"]), 0.95, tolerance = 0.05)
  expect_equal(mean(res[, "cov2"]), 0.95, tolerance = 0.05)
  expect_equal(mean(res[, "cov3"]), 0.95, tolerance = 0.05)

  expect_equal(mean(res[, "bias1"]), 0, tolerance = 0.005)
  expect_equal(mean(res[, "bias2"]), 0, tolerance = 0.005)
  expect_equal(mean(res[, "bias3"]), 0, tolerance = 0.005)

  expect_equal(sd(res[, "est1"]), mean(res[, "se1"]), tolerance = 0.005)
  expect_equal(sd(res[, "est2"]), mean(res[, "se2"]), tolerance = 0.005)
  expect_equal(sd(res[, "est3"]), mean(res[, "se3"]), tolerance = 0.005)
})
