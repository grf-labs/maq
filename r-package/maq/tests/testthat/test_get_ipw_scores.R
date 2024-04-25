test_that("get_ipw_scores works as expected", {
  n <- 5000
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  Y <- -1 +  2 * (W == "B") + 5 * (W == "C") + rnorm(n)
  IPW.scores <- get_ipw_scores(Y, W)
  expect_equal(mean(IPW.scores[, 1]), 2, tolerance = 0.075)
  expect_equal(mean(IPW.scores[, 2]), 5, tolerance = 0.075)

  W.hat <- c(0.2, 0.2, 0.6)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE, prob = W.hat))
  Y <- -1 +  2 * (W == "B") + 5 * (W == "C") + rnorm(n)
  IPW.scores <- get_ipw_scores(Y, W, W.hat = W.hat)
  expect_equal(mean(IPW.scores[, 1]), 2, tolerance = 0.075)
  expect_equal(mean(IPW.scores[, 2]), 5, tolerance = 0.075)

  W <- rbinom(n, 1, 0.5)
  Y <- -1 +  2 * (W == 1) + rnorm(n)
  IPW.scores <- get_ipw_scores(Y, as.factor(W))
  expect_equal(mean(IPW.scores), 2, tolerance = 0.0755)

  W <- rbinom(n, 1, 0.1)
  Y <- -1 +  2 * (W == 1) + rnorm(n)
  IPW.scores <- get_ipw_scores(Y, as.factor(W), W.hat = c(0.9, 0.1))
  expect_equal(mean(IPW.scores), 2, tolerance = 0.1)

  W <- rbinom(n, 1, 0.1)
  Y <- -1 +  2 * (W == 1) + rnorm(n)
  IPW.scores <- get_ipw_scores(Y, as.factor(W), W.hat = c(0.9, 0.1))
  expect_equal(mean(IPW.scores), 2, tolerance = 0.1)
})
