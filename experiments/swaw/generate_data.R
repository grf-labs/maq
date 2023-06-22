generate_data = function(n, p, dgp = c("map")) {
  dgp = match.arg(dgp)
  if (dgp == "map") {
    sigma2 = 4
    A = c(0, 1, 2)
    X = matrix(runif(n * p), n, p)
    R0 = (X[, 5] < 0.6) & (0.35 < X[, 7])
    R2 = (X[, 5]^2 / 0.6^2 + X[, 7]^2 / 0.35^2 < 1) |
      ((X[, 5] - 1)^2 / 0.4^2 + (X[, 7] - 1)^2 / 0.35^2 < 1)
    R1 = !(R0 | R2)
    R = cbind(R0, R1, R2)
    region = which(R, arr.ind = TRUE)
    region = region[order(region[, "row"]), "col"]
    region = region - 1
    
    rewards = sapply(region, function(r) {
      actions <- 0:2
      if (r == 0) {
        return(3 - actions)
      } else if (r == 1) {
        return(2 - abs(actions - 1) / 2)
      } else {
        return(1.5 * (actions - 1))
      }
    })
    rewards = t(rewards)
    colnames(rewards) = 0:2
    realized.a = sample(0:2, n, TRUE)
    reward = rewards[cbind(1:n, realized.a + 1)]
    Y = rnorm(n, mean = reward, sd = sqrt(sigma2))
    
    W = as.factor(realized.a)
    tau = rewards[, -1] - rewards[, 1]
    cost = cbind(X[, 1], 2 * X[, 2])
  }

  list(X = X, Y = Y, W = W, cost = cost, tau = tau)
}

