test_that("maq works as expected", {
  n <- 500
  K <- 5
  reward <- matrix(0.1 + rnorm(n * K), n, K)
  cost <- 0.05 + matrix(runif(n * K), n, K)

  mq <- maq(-abs(reward), cost, reward, R = 150)
  plot(mq)
  summary(mq)
  print(mq)
  plot(mq, add = TRUE, col = 3, ci.args = NULL)
  expect_true(all(predict(mq, 10) == 0))
  expect_equal(average_gain(mq, 10), c(estimate = 0, std.err = 0))

  mq <- maq(reward[, 1], cost[, 1], reward[, 1], R = 100)
  which(predict(mq, 10) > 0)

  # scale invariances
  mq <- maq(reward, cost, reward, seed = 42, R = 100)
  mq.scale <- maq(reward * 1000, cost * 1000, 1000 * reward, seed = 42, R = 100)

  expect_equal(mq[["_path"]]$spend, mq.scale[["_path"]]$spend / 1000, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$gain, mq.scale[["_path"]]$gain / 1000, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$std.err, mq.scale[["_path"]]$std.err / 1000, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$ipath, mq.scale[["_path"]]$ipath)
  expect_equal(mq[["_path"]]$kpath, mq.scale[["_path"]]$kpath)

  # cost scale does not matter
  mq.scalecost <- maq(reward, cost * 1000, reward, seed = 42, R = 100)
  expect_equal(mq[["_path"]]$spend, mq.scalecost[["_path"]]$spend / 1000, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$gain, mq.scalecost[["_path"]]$gain, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$std.err, mq.scalecost[["_path"]]$std.err, tolerance = 1e-12)
  expect_equal(mq[["_path"]]$ipath, mq.scalecost[["_path"]]$ipath)
  expect_equal(mq[["_path"]]$kpath, mq.scalecost[["_path"]]$kpath)

  # incomplete path
  sp <- mq[["_path"]]$spend[100]
  mq.sp <- maq(reward, cost, reward, sp, seed = 42, R = 100)
  len <- length(mq.sp[["_path"]]$spend)
  expect_equal(mq.sp[["_path"]]$spend[len], sp)
  expect_equal(mq[["_path"]]$gain[1:100], mq.sp[["_path"]]$gain)
  expect_equal(mq[["_path"]]$std.err[1:100], mq.sp[["_path"]]$std.err)

  # num.threads don't affect SEs
  mq.1 <- maq(reward, cost, reward, seed = 42, num.threads = 1, R = 100)
  mq.5 <- maq(reward, cost, reward, seed = 42, num.threads = 5, R = 100)
  expect_equal(mq.1[["_path"]]$std.err, mq.5[["_path"]]$std.err)

  # clusters = 1:n
  mq.nocl <- maq(reward, cost, reward, seed = 42, R = 100)
  mq.cl <- maq(reward, cost, reward, seed = 42, clusters = 1:n, R = 100)
  expect_equal(mq.nocl[["_path"]]$std.err, mq.cl[["_path"]]$std.err)
})

test_that("basic invariances hold", {
  n <- 100
  K <- 10
  reward <- matrix(1 + 2 * runif(n * K), n, K)
  cost <- matrix(runif(n * K), n, K)
  reward.eval <- matrix(1 + 2 * runif(n * K), n, K)

  mqini <- maq(reward, cost, reward.eval, seed = 42, R = 200)
  spend <- seq(0.05, 0.7, length.out = 20)

  pi.gain <- c()
  gain <- c()
  for (s in spend) {
    pi.mat <- predict(mqini, spend = s)
    pi.gain <- c(pi.gain, sum(reward.eval * pi.mat) / n)
    gain <-c(gain, average_gain(mqini, spend = s)[[1]])
  }
  expect_equal(pi.gain, gain, tolerance = 1e-12)

  mq.scale <- maq(reward, cost, reward.eval / 10, seed = 42, R = 200)
  expect_equal(mqini[["_path"]]$spend, mq.scale[["_path"]]$spend, tolerance = 1e-12)
  expect_equal(mqini[["_path"]]$gain, mq.scale[["_path"]]$gain * 10, tolerance = 1e-12)
  expect_equal(mqini[["_path"]]$std.err, mq.scale[["_path"]]$std.err * 10, tolerance = 1e-12)
  expect_equal(mqini[["_path"]]$ipath, mq.scale[["_path"]]$ipath)
  expect_equal(mqini[["_path"]]$kpath, mq.scale[["_path"]]$kpath)
})

test_that("sample weighting works as expected", {
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

  mq <- maq(reward, cost, reward, sample.weights = wts, seed = 42, R = 200)
  mq.dupe <- maq(reward.dupe, cost.dupe, reward.dupe, R = 200)

  spends <- c(0.1, 0.25, 0.3, 0.35, 0.4, 0.5)
  est <- lapply(spends, function(s) average_gain(mq, s)[[1]])
  est.dupe <- lapply(spends, function(s) average_gain(mq.dupe, s)[[1]])

  expect_equal(est, est.dupe, tolerance = 1e-12)

  # weight scaling invariance
  mq.scale <- maq(reward, cost, reward, sample.weights = wts * runif(1), seed = 42, R = 200)
  expect_equal(mq[["_path"]], mq.scale[["_path"]], tolerance = 1e-12)
})

test_that("clustering works as expected", {
  n <- 100
  K <- 5
  reward <- matrix(0.1 + rnorm(n * K), n, K)
  cost <- 0.05 + matrix(runif(n * K), n, K)

  rewardc <- rbind(reward, reward, reward, reward, reward)
  costc <- rbind(cost, cost, cost, cost, cost)
  clust <- rep(1:n, 5)

  mq <- maq(reward, cost, reward, R = 200)
  mq.clust <- maq(rewardc, costc, rewardc, clusters = clust, R = 200)

  spends <- c(0.1, 0.25, 0.3, 0.35, 0.4, 0.5)
  est <- lapply(spends, function(s) average_gain(mq, s))
  est.clust <- lapply(spends, function(s) average_gain(mq.clust, s))

  expect_equal(est, est.clust, tolerance = 0.01)
})

test_that("maq with vector cost works as expected", {
  n <- 500
  K <- 5
  reward <- matrix(0.1 + rnorm(n * K), n, K)
  cost <- 0.05 + matrix(runif(n * K), n, K)

  mq1 <- maq(reward[, 1], cost[1], reward[, 1], R = 3)
  mq2 <- maq(reward[, 1], rep(cost[1], n), reward[, 1], R = 3)
  expect_identical(mq1, mq2)

  mq3 <- maq(reward, matrix(cost[1, ], n, K, byrow = TRUE), reward, R = 3)
  mq4 <- maq(reward, cost[1, ], reward, R = 3)
  expect_identical(mq3, mq4)
})

test_that("std errors works as expected", {
  # Get exact coverage of points on the curve
  ntrue <- 100000
  K <- 5
  reward <- matrix(0.1 + rnorm(ntrue * K), ntrue, K)
  cost <- 0.05 + matrix(runif(ntrue * K), ntrue, K)
  mqt <- maq(reward, cost, reward, R = 0)

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
    mq <- maq(reward, cost, reward, R = 150)

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

test_that("paired std errors works as expected", {
  n <- 1000
  K <- 3
  spend <- 0.15

  res <- t(replicate(500, {
    reward <- matrix(0.1 + rnorm(n * K), n, K)
    reward.eval <- matrix(0.1 + rnorm(n * K), n, K)
    cost <- 0.05 + matrix(runif(n * K), n, K)

    mq <- maq(reward, cost, reward.eval, R = 200)
    mq2 <- maq(reward + matrix(0.1 * rnorm(n * K), n, K), cost, reward.eval, R = 200)

    est.diff <- difference_gain(mq, mq2, spend = spend)

    est.1 <- average_gain(mq, spend)
    est.2 <- average_gain(mq2, spend)
    se.naive <- sqrt(est.1[[2]]^2 + est.2[[2]]^2)

    cov <- abs(est.diff[[1]]) / est.diff[[2]] <= 1.96
    cov.naive <- abs(est.diff[[1]]) / se.naive <= 1.96

    c(
      cov = as.numeric(cov),
      cov.naive = as.numeric(cov.naive)
    )
  }))
  cov <- colMeans(res)

  expect_lt(cov[["cov"]], cov[["cov.naive"]])
  expect_equal(cov[["cov"]], 0.95, tolerance = 0.04)
})

test_that("null effect std errors works as expected", {
  n <- 1000
  K <- 10
  s.grid <- seq(0.05, 0.7, length.out = 10)

  res <- t(replicate(500, {
    reward <- matrix(rnorm(n * K), n, K)
    cost <- 0.05 + matrix(runif(n * K), n, K)
    reward.eval <- matrix(rnorm(n * K), n, K)
    mq <- maq(reward, cost, reward.eval, R = 200)

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

test_that("A grf workflow works as expected", {
  n <- 2000
  p <- 5
  X <- matrix(runif(n * p), n, p)
  W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
  Y <- X[, 1] + X[, 2] * (W == "B") + X[, 3] * (W == "C") + rnorm(n)

  spend <- 0.3
  tau.true <- cbind(X[,2], X[,3])
  mq.true <- maq(tau.true, X[, 4:5], tau.true)
  est.true <- average_gain(mq.true, spend = spend)[[1]]

  train <- sample(1:n, n/2)
  eval <- -train
  cost.hat <- X[eval, 4:5]
  tau.forest <- grf::multi_arm_causal_forest(X[train, ], Y[train], W[train])
  tau.hat <- predict(tau.forest, X[eval, ])$predictions[,,]

  # A MAQ evaluated via DR scores
  eval.forest <- grf::multi_arm_causal_forest(X[eval, ], Y[eval], W[eval])
  DR.scores <- grf::get_scores(eval.forest)[,,]
  mq.aipw <- maq(tau.hat, cost.hat, DR.scores, R = 200)
  est.aipw <- average_gain(mq.aipw, spend = spend)

  # A MAQ evaluated via IPW
  W.hat.true <- rep(1/3, 3)
  observed.W <- match(W, levels(W))
  Y.k.mat <- matrix(0, length(W), nlevels(W))
  Y.k.mat[cbind(seq_along(observed.W), observed.W)] <- Y
  Y.k.ipw <- sweep(Y.k.mat, 2, W.hat.true, "/")
  Y.k.ipw.eval <- Y.k.ipw[eval, -1] - Y.k.ipw[eval, 1]

  mq.ipw <- maq(tau.hat, cost.hat, Y.k.ipw.eval, R = 200)
  est.ipw <- average_gain(mq.ipw, spend = spend)

  expect_equal(est.aipw[[1]], est.ipw[[1]], tolerance = 0.15)
  expect_equal(est.aipw[[1]], est.true, tolerance = 3 * est.aipw[[2]])
  expect_equal(est.ipw[[1]], est.true, tolerance = 3 * est.ipw[[2]])
})

test_that("tie handling works as expected", {
  n <- 500
  K <- 3
  reward <- matrix(sample(c(1, 10), n * K, TRUE), n, K)
  cost <- matrix(sample(c(1, 2), n * K, TRUE), n, K)

  mq1 <- maq(reward, cost, reward, R = 0)
  mq2 <- maq(reward, cost, reward, tie.breaker = rev(1:n), R = 0)

  expect_true(all(mq1[["_path"]]$ipath[1:20] < mq2[["_path"]]$ipath[1:20]))
})

test_that("avg maq works as expected", {
  # Same point estimates as "duplicated" maq-approach
  n <- 1000
  K <- 15
  reward <- matrix(0.1 + rnorm(n * K), n, K)
  reward.eval <- matrix(0.1 + rnorm(n * K), n, K)
  cost <- 1:K

  mqr <- maq(reward, cost, reward.eval, target.with.covariates = FALSE)
  mq <- maq(matrix(colMeans(reward), n, K, byrow = TRUE),
            matrix(cost, n, K, byrow = TRUE),
            matrix(colMeans(reward.eval), n, K, byrow = TRUE))
  expect_equal(mqr[["_path"]]$spend, mq[["_path"]]$spend)
  expect_equal(mqr[["_path"]]$gain, mq[["_path"]]$gain)

  # < 0
  reward <- matrix(-10 + rnorm(n * K), n, K)
  mqr <- maq(reward, cost, reward.eval, target.with.covariates = FALSE)
  mq <- maq(matrix(colMeans(reward), n, K, byrow = TRUE),
            matrix(cost, n, K, byrow = TRUE),
            matrix(colMeans(reward.eval), n, K, byrow = TRUE))
  expect_equal(mqr[["_path"]]$spend, mq[["_path"]]$spend)
  expect_equal(mqr[["_path"]]$gain, mq[["_path"]]$gain)

  # std.errors works as expected
  ntrue <- 100000
  K <- 5
  rewardt <- cbind(runif(ntrue), matrix(-0.1 + rnorm(ntrue * K), ntrue, K))
  cost <- c(0.1, rep(1, K))
  mqt <- maq(rewardt, cost, rewardt, target.with.covariates = FALSE)

  spend1 <- 0.02
  true1 <- average_gain(mqt, spend1)[[1]]
  spend2 <- 0.1
  true2 <- average_gain(mqt, spend2)[[1]]

  n <- 500
  res <- t(replicate(250, {
    reward <- rewardt[sample(ntrue, n), ]
    mq <- maq(reward, cost, reward, target.with.covariates = FALSE, R = 200)

    est1 <- average_gain(mq, spend1)
    est2 <- average_gain(mq, spend2)

    c(
      est1 = est1[[1]],
      se1 = est1[[2]],
      bias1 = est1[[1]] - true1,
      cov1 = abs(est1[[1]] - true1) / est1[[2]] <= 1.96,

      est2 = est2[[1]],
      se2 = est2[[2]],
      bias2 = est2[[1]] - true2,
      cov2 = abs(est2[[1]] - true2) / est2[[2]] <= 1.96
    )
  }))

  expect_equal(mean(res[, "cov1"]), 0.95, tolerance = 0.04)
  expect_equal(mean(res[, "cov2"]), 0.95, tolerance = 0.04)

  # std.errors works as expected (larger avg hull)
  rewardt <- cbind(runif(ntrue), 2 * runif(ntrue), 3 * runif(ntrue), -1 * runif(ntrue))
  cost <- c(0.001, 0.02, 0.06, 1)
  mqt <- maq(rewardt, cost, rewardt, target.with.covariates = FALSE)

  spend1 <- 0.002
  true1 <- average_gain(mqt, spend1)[[1]]
  spend2 <- 0.02
  true2 <- average_gain(mqt, spend2)[[1]]
  spend3 <- 0.05
  true3 <- average_gain(mqt, spend3)[[1]]

  n <- 500
  res <- t(replicate(250, {
    reward <- rewardt[sample(ntrue, n), ]
    mq <- maq(reward, cost, reward, target.with.covariates = FALSE, R = 200)

    est1 <- average_gain(mq, spend1)
    est2 <- average_gain(mq, spend2)
    est3 <- average_gain(mq, spend3)

    c(
      est1 = est1[[1]],
      se1 = est1[[2]],
      bias1 = est1[[1]] - true1,
      cov1 = abs(est1[[1]] - true1) / est1[[2]] <= 1.96,

      est2 = est2[[1]],
      se2 = est2[[2]],
      bias2 = est2[[1]] - true2,
      cov2 = abs(est2[[1]] - true2) / est2[[2]] <= 1.96,

      est3 = est3[[1]],
      se3 = est3[[2]],
      bias3 = est3[[1]] - true3,
      cov3 = abs(est3[[1]] - true3) / est3[[2]] <= 1.96
    )
  }))

  expect_equal(mean(res[, "cov1"]), 0.95, tolerance = 0.04)
  expect_equal(mean(res[, "cov1"]), 0.95, tolerance = 0.04)
  expect_equal(mean(res[, "cov3"]), 0.95, tolerance = 0.04)
})

test_that("predict type works as expected", {
  n <- 500
  K <- 10
  reward <- matrix(rnorm(n * K), n, K)
  cost <- matrix(runif(n * K), n, K)
  DR.scores <- reward + rnorm(n)
  mq <- maq(reward, cost, DR.scores)

  pi.mat <- predict(mq, 0.1)
  pi.vec <- predict(mq, 0.1, type = "vector")

  expect_lte(
    sum(cost[cbind(1:n, pi.vec)]) / n,
    sum(cost * pi.mat) / n
  )

  data.path <- summary(mq)
  some.value <- sample(nrow(data.path), 1)
  spend <- data.path$spend[some.value]

  pi.mat.int <- predict(mq, spend)
  pi.vec.int <- predict(mq, spend, type = "vector")

  expect_equal(
    sum(cost[cbind(1:n, pi.vec.int)]) / n,
    sum(cost * pi.mat.int) / n
  )
  expect_equal(predict(mq, -10, type = "vector"), rep(0, n))
})
