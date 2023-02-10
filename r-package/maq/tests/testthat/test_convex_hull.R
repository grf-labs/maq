# TODO paste cvx hull function here

test_that("convex hull works as expected", {
  # simple test cases
  cost <- rbind(
    c(5, 0.5, 3, 1, 4),
    c(5, 0.5, 3, 1, 4),
    c(5, 0.5, 3, 1, 4),
    c(5, 0.5, 3, 1, 4),
    c(5, 0.5, 3, 1, 4),
    c(2, 2, 3, 1, 4),
    c(4, 2, 3, 1, 4.2)
  )
  reward <- rbind(
    c(3, 2, 4, 2, 3.9),
    c(3, 2, 3.8, 2, 3.9),
    c(3, 0.5, 3.8, 2, 3.9),
    c(3.91, 2, 3.8, 2, 3.9),
    c(3.9, 2, 3.8, 2, 3.9),
    c(3, 2, 4, 2, 3.9),
    c(4, 2, 4, 2.3, 4)
  )
  R.expected <- list(
    c(2, 3),
    c(2, 3, 5),
    c(4, 3, 5),
    c(2, 3, 5, 1),
    c(2, 3, 5),
    c(4, 1, 3),
    c(4, 3)
  )
  R <- convex_hull(reward, cost)
  expect_equal(R, R.expected)

  # /w negative rewards
  cost <- rbind(
    c(4, 2, 3, 1, 4.2),
    c(4, 2, 3, 1, 4.2),
    c(4, 2, 3, 1, 4.2),
    c(4, 2, 3, 1, 4.2),
    c(4, 2, 3, 1, 4.2),
    c(4, 2, 3, 1, 4.2)
  )
  reward <- rbind(
    c(4, 2, -2, 2.3, 4),
    c(4, 2, -2, -1, 4),
    c(-4, -2, 2, -1, -4),
    c(4, 1.9, 2, -2, 4),
    c(-4, -1.9, -2, -2, 4),
    c(-4, -1.2, -2, -1, -1)
  )
  R.neg.expected <- list(
    c(4, 1),
    c(2, 1),
    c(3),
    c(1),
    c(5),
    numeric(0)
  )
  R.neg <- convex_hull(reward, cost)
  expect_equal(R.neg, R.neg.expected)

  # various
  expect_equal(
    convex_hull(reward = rbind(c(-1, -1)), cost = rbind(c(1, 1))),
    list(numeric(0))
  )
  expect_equal(
    convex_hull(reward = rbind(c(-1, -10)), cost = rbind(c(1, 0.5))),
    list(numeric(0))
  )
  expect_equal(
    convex_hull(reward = rbind(c(0, 0)), cost = rbind(c(1, 0.5))),
    list(numeric(0))
  )
  expect_equal(
    convex_hull(reward = rbind(c(0, 10)), cost = rbind(c(1, 0.5))),
    list(2)
  )
  expect_equal(
    convex_hull(reward = rbind(c(10, 10.1)), cost = rbind(c(1, 1))),
    list(2)
  )
  expect_equal(
    convex_hull(reward = rbind(c(-10, -10.1)), cost = rbind(c(0, 0))),
    list(numeric(0))
  )
  expect_equal(
    convex_hull(reward = rbind(c(10, 10, 10, 10)), cost = rbind(c(1, 0.5, 10, 0.1))),
    list(4)
  )
  expect_equal(
    convex_hull(reward = rbind(c(-13, 10, 12)), cost = rbind(c(0.5, 1, 1))),
    list(3)
  )
  expect_equal(
    convex_hull(reward = rbind(c(-13, 12, 10)), cost = rbind(c(0.5, 1, 1))),
    list(2)
  )
  expect_equal(
    convex_hull(reward = rbind(c(-13, 12, 10, 19, -9)), cost = rbind(c(0.5, 1, 1, 1, 1))),
    list(4)
  )
  expect_equal(
    convex_hull(reward = rbind(c(-13, -13, -13)), cost = rbind(c(1, 1, 1))),
    list(numeric(0))
  )
  expect_equal(
    convex_hull(reward = rbind(c(-13, -13, 13)), cost = rbind(c(1, 1, 1))),
    list(3)
  )

  # stable wrt cost sort order
  expect_equal(
    convex_hull(reward = rbind(c(-13, 13, 13)), cost = rbind(c(1, 1, 1))),
    list(2)
  )
  expect_equal(
    convex_hull(reward = rbind(c(13, -13, 13)), cost = rbind(c(1, 1, 1))),
    list(1)
  )
  expect_equal(
    convex_hull(reward = rbind(c(13, 13, 13)), cost = rbind(c(1, 1, 1))),
    list(1)
  )
  expect_equal(
    convex_hull(reward = rbind(c(-13, 13, 14)), cost = rbind(c(1, 1, 1))),
    list(3)
  )
  expect_equal(
    convex_hull(reward = rbind(c(13, 13, 13, 14, 14)), cost = rbind(c(1, 1, 1, 2, 2))),
    list(c(1, 4))
  )

  # basic invariances
  n <- 500
  K <- 35
  reward <- matrix(1 + rnorm(n * K), n, K)
  cost <- 0.05 + matrix(runif(n * K), n, K)
  expect_equal(
    convex_hull(reward, cost),
    convex_hull(reward/1000, cost)
  )
  expect_equal(
    convex_hull(reward, cost),
    convex_hull(reward/1000, cost/100)
  )
  expect_equal(
    unlist(convex_hull(-abs(reward), cost)),
    numeric(0)
  )

})
