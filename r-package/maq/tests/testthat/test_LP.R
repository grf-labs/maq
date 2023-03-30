# Solves the linear knapsack-type (MCKP) problem with a generic LP solver for a given spend.
tryCatch(
  {
    attachNamespace("lpSolve")
  },
  error = function(e) {
    install.packages("lpSolve", repos = "http://cran.us.r-project.org")
    attachNamespace("lpSolve")
  }
)
lp_solver = function(reward, cost, budget) {
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

  mqini <- maq(reward, cost, budget, reward)
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

  mqini <- maq(reward, cost, budget, reward)

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

  mqini <- maq(reward, cost, budget, reward)
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

  mqini <- maq(reward, cost, budget, reward)
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

test_that("non-unique solution works as expected", {
  budget <- 100
  n <- 100
  K <- 5
  reward <- matrix(sample(c(1, 10), n * K, TRUE), n, K)
  cost <- matrix(sample(c(1, 2), n * K, TRUE), n, K)
  spend <- 1

  mq1 <- maq(reward, cost, budget, reward, R = 0)
  mq2 <- maq(reward, cost, budget, reward, R = 0, tie.breaker = sample(1:n))
  lp <- lp_solver(reward, cost, spend)
  lp.reward <- sum(lp$alloc.mat * reward) / n

  expect_equal(average_gain(mq1, spend = spend)[[1]], lp.reward, tolerance = 1e-10)
  expect_equal(average_gain(mq2, spend = spend)[[1]], lp.reward, tolerance = 1e-10)
})

test_that("SEs capture LP re-solved", {
  budget <- 100
  n <- 100
  K <- 3
  reward <- matrix(rnorm(n * K), n, K)
  cost <- 0.05 + matrix(runif(n * K), n, K)
  R <- 500
  mq <- maq(reward, cost, budget, reward, R = R)

  # pick an arbitrary spend point except the initial grid point
  sp <- sample(mq[["_path"]]$spend[-(1:3)], 1)
  mq.se <- average_gain(mq, spend = sp)[[2]]
  # this SE should correspond to

  clusters <- 1:n
  samples.by.cluster <- split(seq_along(clusters), clusters)
  n.bs <- floor(length(samples.by.cluster) / 2)
  index.list <- replicate(R, unlist(samples.by.cluster[sample.int(n, n.bs, replace = FALSE)], use.names = FALSE), simplify = FALSE)

  resbs <- lapply(index.list, function(ix) {
    lp_solver(reward[ix, ], cost[ix, ], sp)$gain
  })
  expect_equal(mq.se, sd(unlist(resbs)), tolerance = 0.005)

  # same, with clusters
  clusters <- rep(1:10, 10)
  mq.cl <- maq(reward, cost, 100, reward, R = R, clusters = clusters)
  mq.se.cl <- average_gain(mq.cl, sp)[[2]]

  samples.by.cluster <- split(seq_along(clusters), clusters)
  n <- length(samples.by.cluster)
  n.bs <- floor(length(samples.by.cluster) / 2)
  index.list <- replicate(R, unlist(samples.by.cluster[sample.int(n, n.bs, replace = FALSE)], use.names = FALSE), simplify = FALSE)

  resbs <- lapply(index.list, function(ix) {
    lp_solver(reward[ix, ], cost[ix, ], sp)$gain
  })
  expect_equal(mq.se.cl, sd(unlist(resbs)), tolerance = 0.005)
})
