rm(list = ls())
set.seed(123)
source("generate_data.R")
library(maq)
library(grf)

n = 5000
p = 10
data.train = generate_data(n, p)
data.eval = generate_data(n, p)

cf.train = multi_arm_causal_forest(
  data.train$X,
  data.train$Y,
  data.train$W
)
tau.hat.eval = predict(cf.train, data.eval$X, drop = TRUE)$predictions

cf.eval = multi_arm_causal_forest(
  data.eval$X,
  data.eval$Y,
  data.eval$W
)
dr.eval = get_scores(cf.eval, drop = TRUE)
cost.eval = data.eval$cost

mq = maq(tau.hat.eval, cost.eval, dr.eval, R = 200)
mq.1 = maq(tau.hat.eval[, 1], cost.eval[, 1], dr.eval[, 1], R = 200)
mq.2 = maq(tau.hat.eval[, 2], cost.eval[, 2], dr.eval[, 2], R = 200)
mq.r = maq(tau.hat.eval, cost.eval, dr.eval, target.with.covariates = FALSE, R = 200)

pdf("figure_1.pdf")
plot(mq, lwd = 2, ci.args = NULL)
plot(mq.1, add = TRUE, lty = 2, col = 2, lwd = 2, ci.args = NULL)
plot(mq.2, add = TRUE, lty = 2, col = 3, lwd = 2, ci.args = NULL)
plot(mq.r, add = TRUE, lty = 1, col = 4, lwd = 2, ci.args = NULL)
legend("topleft", c("Qini (multi-armed)", "Qini (arm 1)", "Qini (arm 2)", "Without targeting"), col = 1:4, lty = c(1,2,2,1), lwd = 2)
dev.off()

# std.err
spend = 0.5
options(digits=2)
average_gain(mq.r, spend)
average_gain(mq.1, spend)
average_gain(mq, spend)
