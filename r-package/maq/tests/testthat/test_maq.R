test_that("maq works as expected", {
  budget <- 1000
  n <- 500
  K <- 5
  reward <- matrix(0.1 + rnorm(n * K), n, K)
  cost <- 0.05 + matrix(runif(n * K), n, K)

  mq <- maq(-abs(reward), cost, budget, reward, R = 150)
  plot(mq)
  summary(mq)
  print(mq)
  expect_true(all(predict(mq, 10) == 0))
  expect_equal(average_gain(mq, 10), c(estimate = 0, std.err = 0))

  mq <- maq(reward[, 1], cost[, 1], budget, reward[, 1], R = 150)
  Matrix::which(predict(mq, 10) > 0)

  # scale invariances
  mq <- maq(reward, cost, budget, reward, seed = 42)
  mq.scale <- maq(reward * 1000, cost * 1000, budget, 1000 * reward, seed = 42)

  expect_equal(mq[["_path"]]$spend, mq.scale[["_path"]]$spend / 1000, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$gain, mq.scale[["_path"]]$gain / 1000, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$std.err, mq.scale[["_path"]]$std.err / 1000, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$ipath, mq.scale[["_path"]]$ipath)
  expect_equal(mq[["_path"]]$kpath, mq.scale[["_path"]]$kpath)

  # cost scale does not matter
  mq.scalecost <- maq(reward, cost * 1000, budget, reward, seed = 42)
  expect_equal(mq[["_path"]]$spend, mq.scalecost[["_path"]]$spend / 1000, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$gain, mq.scalecost[["_path"]]$gain, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$std.err, mq.scalecost[["_path"]]$std.err, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$ipath, mq.scalecost[["_path"]]$ipath)
  expect_equal(mq[["_path"]]$kpath, mq.scalecost[["_path"]]$kpath)
})

test_that("basic invariances hold", {
  budget <- 25
  n <- 100
  K <- 10
  reward <- matrix(1 + 2 * runif(n * K), n, K)
  cost <- matrix(runif(n * K), n, K)
  reward.eval <- matrix(1 + 2 * runif(n * K), n, K)

  mqini <- maq(reward, cost, budget, reward.eval, seed = 42)
  spend <- seq(0.05, 0.7, length.out = 20)

  pi.gain <- c()
  gain <- c()
  for (s in spend) {
    pi.mat <- predict(mqini, spend = s)
    pi.gain <- c(pi.gain, sum(reward.eval * pi.mat) / n)
    gain <-c(gain, average_gain(mqini, spend = s)[[1]])
  }
  expect_equal(pi.gain, gain, tolerance = 1e-12)

  mq.scale <- maq(reward, cost, budget, reward.eval / 10, seed = 42)
  expect_equal(mqini[["_path"]]$spend, mq.scale[["_path"]]$spend, tolerance = 1e-12)
  expect_equal(mqini[["_path"]]$gain, mq.scale[["_path"]]$gain * 10, tolerance = 1e-12)
  expect_equal(mqini[["_path"]]$std.err, mq.scale[["_path"]]$std.err * 10, tolerance = 1e-12)
  expect_equal(mqini[["_path"]]$ipath, mq.scale[["_path"]]$ipath)
  expect_equal(mqini[["_path"]]$kpath, mq.scale[["_path"]]$kpath)
})

test_that("sample weighting works as expected", {
  budget <- 1000
  n <- 200
  K <- 5
  reward <- matrix(0.1 + rnorm(n * K), n, K)
  cost <- 0.05 + matrix(runif(n * K), n, K)

  # giving weight 2 exactly same sample as duplicating
  dupe <- sample(1:n, 100)
  reward.dupe <- rbind(reward, reward[dupe, ])
  cost.dupe <- rbind(cost, cost[dupe, ])
  wts <- rep(1, n)
  wts[dupe] <- 2

  mq <- maq(reward, cost, budget, reward, sample.weights = wts, seed = 42)
  mq.dupe <- maq(reward.dupe, cost.dupe, budget, reward.dupe)

  spends <- c(0.1, 0.25, 0.3, 0.35, 0.4, 0.5)
  est <- lapply(spends, function(s) average_gain(mq, s)[[1]])
  est.dupe <- lapply(spends, function(s) average_gain(mq.dupe, s)[[1]])

  expect_equal(est, est.dupe, tolerance = 1e-12)

  # weight scaling invariance
  mq.scale <- maq(reward, cost, budget, reward, sample.weights = wts * runif(1), seed = 42)
  expect_equal(mq[["_path"]], mq.scale[["_path"]], tolerance = 1e-12)
})

test_that("clustering works as expected", {
  budget <- 1000
  n <- 100
  K <- 5
  reward <- matrix(0.1 + rnorm(n * K), n, K)
  cost <- 0.05 + matrix(runif(n * K), n, K)

  rewardc <- rbind(reward, reward, reward, reward, reward)
  costc <- rbind(cost, cost, cost, cost, cost)
  clust <- rep(1:n, 5)

  mq <- maq(reward, cost, budget, reward)
  mq.clust <- maq(rewardc, costc, budget, rewardc, clusters = clust)

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
  mqt <- maq(reward, cost, budget, reward, R = 0)

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
    mq <- maq(reward, cost, budget, reward, R = 150)

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

test_that("null effect std errors works as expected", {
  budget <- 100
  n <- 1000
  K <- 10
  s.grid <- seq(0.05, 0.7, length.out = 10)

  res <- t(replicate(500, {
    reward <- matrix(rnorm(n * K), n, K)
    cost <- 0.05 + matrix(runif(n * K), n, K)
    reward.eval <- matrix(rnorm(n * K), n, K)
    mq <- maq(reward, cost, budget, reward.eval)

    est <- lapply(s.grid, function(s) average_gain(mq, s))
    df.est <- do.call(rbind, est)

    z.stat <- abs(df.est[, "estimate"] - 0) / df.est[, "std.err"]
    coverage <- z.stat <= 1.96
    coverage[is.na(coverage)] <- 1 # all cates zero?

    c(
      est = df.est[, "estimate"],
      se = df.est[, "std.err"],
      cov = coverage
    )
  }))
  iest <- 1:length(s.grid)
  ise <- iest + 10
  icov <- -c(iest, ise)

  expect_equal(
    unname(colMeans(res[, iest])),
    rep(0, length(s.grid)),
    tolerance = 0.01
  )
  expect_equal(
    unname(colMeans(res[, ise])),
    unname(apply(res[, iest], 2, sd)),
    tolerance = 0.01
  )
  expect_equal(
    unname(colMeans(res[, icov])),
    rep(0.95, length(s.grid)),
    tolerance = 0.05
  )
})

test_that("tie handling works as expected", {
  budget <- 100
  n <- 500
  K <- 3
  reward <- matrix(sample(c(1, 10), n * K, TRUE), n, K)
  cost <- matrix(sample(c(1, 2), n * K, TRUE), n, K)

  mq1 <- maq(reward, cost, budget, reward, R = 0)
  mq2 <- maq(reward, cost, budget, reward, tie.breaker = rev(1:n), R = 0)

  expect_true(all(mq1[["_path"]]$ipath[1:20] < mq2[["_path"]]$ipath[1:20]))
})
