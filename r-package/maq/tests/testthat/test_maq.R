test_that("std errors work as expected", {
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
