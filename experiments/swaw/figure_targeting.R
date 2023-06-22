rm(list = ls())
set.seed(42)
library(maq)
library(grf)

noise_sd = 1.5
base_mu = 1
cost_covariate = function(x, k) {
  if (k == 1) {
    return (2.5*(x < 0.3) + 2*(x >= 0.3)*(x <= 0.7) + 3.5*(x > 0.7))
  } else if (k == 2) {
    return (3*(x < 0.3) + 2.5*(x >= 0.3)*(x <= 0.7) + 3.0*(x > 0.7))
  } else {
    return (2*(x < 0.3) + 3*(x >= 0.3)*(x <= 0.7) + 4*(x > 0.7))
  }
}

outcome_covariate = function(x, k) {
  if (k == 0) {
    return (base_mu)
  }else if (k == 1) {
    return (3.35*(x < 0.3) + 3.42*(x >= 0.3)*(x <= 0.7) + 4.5*(x > 0.7))
  } else if (k ==2) {
    return (8.2*(x < 0.2) + 1.6*(x >= 0.2)*(x <= 0.7) + 2.1*(x > 0.7))
  } else {
    return (2.5*(x < 0.3) + 2.0*(x >= 0.3)*(x <= 0.6) + 6*(x > 0.6))
  }
}

generate_data = function(n) {
  W = sample(0:3, size=n, replace = TRUE)
  X = as.matrix(runif(n, 0, 1))
  mu = cbind(sapply(1:n, function(i) outcome_covariate(X[i], 1)),
             sapply(1:n, function(i) outcome_covariate(X[i], 2)),
             sapply(1:n, function(i) outcome_covariate(X[i], 3)))
  tau = mu - base_mu
  Y = numeric(n)
  for (i in 1:n) {
    if (W[i] != 0) {
      Y[i] = mu[i, W[i]]
    } else {
      Y[i] = base_mu
    }
  }
  Y = Y + rnorm(n, 0, noise_sd)

  cost_oracle = cbind(sapply(1:n, function(i) cost_covariate(X[i], 1)),
                      sapply(1:n, function(i) cost_covariate(X[i], 2)),
                      sapply(1:n, function(i) cost_covariate(X[i], 3)))

  cost = numeric(n)
  for (i in 1:n) {
    if (W[i] != 0) {
     cost[i] = cost_oracle[i,W[i]]+rnorm(1,0,noise_sd)
    }
  }
  W = as.factor(W)
  list(X = X, Y = Y, W = W, cost = cost, cost.oracle = cost_oracle, tau = tau)
}

max.budget = 10
n_eval = 1500
data.eval = generate_data(n_eval)
cf.eval = multi_arm_causal_forest(
   data.eval$X,
   data.eval$Y,
   data.eval$W
)
dr.eval = get_scores(cf.eval, drop = TRUE)

data.eval$tau = matrix(rnorm(n_eval*3,0,1.2),n_eval, 3) + data.eval$tau
data.eval$cost.oracle = pmax(matrix(0.001, n_eval,3), matrix(rnorm(n_eval*3,0,1.2),n_eval, 3) + data.eval$cost.oracle)

mq = maq(data.eval$tau, data.eval$cost.oracle, max.budget, dr.eval)
mq.1 = maq(data.eval$tau[, 1], data.eval$cost.oracle[,1], max.budget, dr.eval[, 1])
mq.2 = maq(data.eval$tau[, 2], data.eval$cost.oracle[,2], max.budget, dr.eval[, 2])
mq.3 = maq(data.eval$tau[, 3], data.eval$cost.oracle[,3], max.budget, dr.eval[, 3])
mq.12 = maq(data.eval$tau[, c(1,2)], data.eval$cost.oracle[,c(1,2)], max.budget, dr.eval[, c(1,2)])
mq.13 = maq(data.eval$tau[, c(1,3)], data.eval$cost.oracle[,c(1,3)], max.budget, dr.eval[, c(1,3)])
mq.23 = maq(data.eval$tau[, c(2,3)], data.eval$cost.oracle[,c(2,3)], max.budget, dr.eval[, c(2,3)])
mq.random = maq(data.eval$tau, data.eval$cost.oracle, max.budget, dr.eval, target.with.covariates = FALSE)
mq.random1 = maq(data.eval$tau[, 1], data.eval$cost.oracle[, 1], max.budget, dr.eval[,1], target.with.covariates = FALSE)
mq.random2 = maq(data.eval$tau[, 2], data.eval$cost.oracle[, 2], max.budget, dr.eval[,2], target.with.covariates = FALSE)
mq.random3 = maq(data.eval$tau[, 3], data.eval$cost.oracle[, 3], max.budget, dr.eval[,3], target.with.covariates = FALSE)


# Panel A
pdf("figure_targeting_a.pdf")
plot(mq, col = "red", ci.args = NULL, lwd = 2,cex.sub=1.6, cex.axis=1.4, cex.lab=1.4)
plot(mq.random, add = TRUE, col = "coral", lwd = 2)
plot(mq.random1, add = TRUE, col = "blue", lwd = 2)
plot(mq.random2, add = TRUE, col = "purple", lwd = 2)
plot(mq.random3, add = TRUE, col = "green", lwd = 2)
legend("topleft",
       c("All arms (targeting)", "All arms (no targeting)", "Arm 1 (no targeting)", "Arm 2 (no targeting)", "Arm 3 (no targeting"),
       col = c("red", "coral", "blue", "purple", "green"), lty = 1, lwd = 2)
dev.off()

# Panel B
pdf("figure_targeting_b.pdf")
plot(mq, col = "red", lwd = 2,cex.sub=1.6, cex.axis=1.4, cex.lab=1.4)
plot(mq.1, add = TRUE, col = "blue", lwd = 2)
plot(mq.2, add = TRUE, col = "purple", lwd = 2)
plot(mq.3, add = TRUE, col = "green", lwd = 2)
plot(mq.12, add = TRUE, col = "brown", lwd = 2)
plot(mq.13, add = TRUE, col = "black", lwd = 2)
plot(mq.23, add = TRUE, col = "dark gray", lwd = 2)
legend("topleft",
       c("All arms", "Arm 1", "Arm 2", "Arm 3", "Arm 1,2", "Arm 1,3", "Arm 2,3"),
       col = c("red", "blue", "purple", "green", "brown", "black", " dark gray"), lty = 1, lwd = 2)
dev.off()

