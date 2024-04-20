set.seed(123)
rm(list = ls())
library(maq)
library(grf)
source("generate_data.R")

data.train = generate_data(n = 5000, p = 10)
data.eval = generate_data(n = 4000, p = 10)

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

mq = maq(tau.hat.eval, cost.eval, dr.eval)
mq.1 = maq(tau.hat.eval[, 1], cost.eval[, 1], dr.eval[, 1], R = 200)
mq.2 = maq(tau.hat.eval[, 2], cost.eval[, 2], dr.eval[, 2], R = 200)
mq.r = maq(tau.hat.eval, cost.eval, dr.eval, target.with.covariates = FALSE, R = 200)

pdf("figure_1.pdf", pointsize = 16)
plot(mq, lwd = 2, ci.args = NULL)
plot(mq.1, add = TRUE, lty = 2, col = 2, lwd = 3, ci.args = NULL)
plot(mq.2, add = TRUE, lty = 3, col = 4, lwd = 4, ci.args = NULL)
plot(mq.r, add = TRUE, lty = 4, col = 8, lwd = 4, ci.args = NULL)
legend("topleft",
       c("Multi-armed", "Only arm 1", "Only arm 2", "No targeting"),
       col = c(1, 2, 4, 8), lty = 1:4, lwd = 2, bty = "n")
dev.off()

# Estimates
options(digits=1)
spend = 0.2
average_gain(mq.1, spend)
average_gain(mq.2, spend)
average_gain(mq, spend)
