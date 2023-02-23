# TODO: paste lp_solver func in here

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
    # spend <- spend * n
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
  # spend <- spend * n
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
