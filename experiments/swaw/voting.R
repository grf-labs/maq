rm(list = ls())
set.seed(123)
library(maq)
library(grf)

# Data file "D001F01" available to download from https://isps.yale.edu/research/data/d001 and contains "GerberGreenLarimer_APSR_2008_social_pressure.csv"
raw = read.csv("GerberGreenLarimer_APSR_2008_social_pressure.csv")
raw$treatment = trimws(raw$treatment)
# length(unique(raw$hh_id))

cluster = as.integer(as.factor(raw$hh_id))
Y = as.integer(raw$voted == "Yes")
W = factor(raw$treatment, levels = c("Control", "Civic Duty", "Hawthorne", "Self", "Neighbors"))
X = cbind(
  year.of.birth = raw$yob,
  sex = as.integer(raw$sex == "male"),
  household.size = raw$hh_size,
  voted.general2000 = as.integer(raw$g2000 == "yes"),
  voted.general2002 = as.integer(raw$g2002 == "yes"),
  voted.general2004 = as.integer(raw$g2004 == "yes"),
  voted.primary2000 = as.integer(raw$p2000 == "yes"),
  voted.primary2002 = as.integer(raw$p2002 == "yes"),
  voted.primary2004 = as.integer(raw$p2004 == "Yes")
)

# The true propensities
W.hat = c(10/18, 2/18, 2/18, 2/18, 2/18)
# summary(W) / length(W)

# Draw a random sample of households for train/evaluation (50/50)
samples.by.hh = split(seq_along(cluster), cluster)
num.hh = length(samples.by.hh)
train = unlist(samples.by.hh[sample(1:num.hh, num.hh / 2)], use.names = FALSE)
eval = -train

cf.train = multi_arm_causal_forest(X[train, ], Y[train], W[train], W.hat = W.hat, clusters = cluster[train])
tau.hat.eval = predict(cf.train, X[eval, ], drop = TRUE)$predictions
average_treatment_effect(cf.train)
head(tau.hat.eval)

# Evaluate via IPW
observed.W = match(W, levels(W))
Y.k.mat = matrix(0, length(W), nlevels(W))
Y.k.mat[cbind(seq_along(observed.W), observed.W)] = Y
Y.k.ipw = sweep(Y.k.mat, 2, W.hat, "/")
Y.k.ipw.eval = Y.k.ipw[eval, -1] - Y.k.ipw[eval, 1]

## save/load R-session with fit forests
# save.image("voting.RData", compress = TRUE)
# load("voting.RData")
##

# Fit cost curves
cost = c(1, 15, 30, 45)
mq = maq(tau.hat.eval, cost, Y.k.ipw.eval, clusters = cluster[eval], R = 200)
mq.civic = maq(tau.hat.eval[, 1], cost[1], Y.k.ipw.eval[, 1], clusters = cluster[eval], R = 200)
mq.hawthorne = maq(tau.hat.eval[, 2], cost[2], Y.k.ipw.eval[, 2], clusters = cluster[eval], R = 200)
mq.self = maq(tau.hat.eval[, 3], cost[3], Y.k.ipw.eval[, 3], clusters = cluster[eval], R = 200)
mq.neighbors = maq(tau.hat.eval[, 4], cost[4], Y.k.ipw.eval[, 4], clusters = cluster[eval], R = 200)
# Non-targeting baseline
mq.avg = maq(tau.hat.eval, cost, Y.k.ipw.eval, clusters = cluster[eval], R = 200, target.with.covariates = FALSE)


spend = 5
# Plot cost curves
pdf("voting_figure_a.pdf")
plot(mq, xlab = "Intrusion cost", ylab = "Increase in probability of voting", ci.args = NULL,  ylim = c(0, 0.06), lwd = 2)
plot(mq.civic, add = TRUE, col = 2, ci.args = NULL, lwd = 2)
plot(mq.hawthorne, add = TRUE, col = 3, ci.args = NULL, lwd = 2)
plot(mq.self, add = TRUE, col = 4, ci.args = NULL, lwd = 2)
plot(mq.neighbors, add = TRUE, col = 7, ci.args = NULL, lwd = 2)
abline(v = spend, lty = 2)
legend("topleft", c("All", "Civic", "Hawthorne", "Self", "Neighbors"), col = c(1,2,3,4,7), lty = 1, bg = "white", lwd = 2)
dev.off()

pdf("voting_figure_b.pdf")
plot(mq, xlab = "Intrusion cost", ylab = "Increase in probability of voting",  ylim = c(0, 0.06), lwd = 2)
plot(mq.avg, add = TRUE, col = 6, lwd = 2)
abline(v = spend, lty = 2)
legend("topleft", c("All arms (with targeting)", "All arms (without targeting)", "95% CI"), col = c(1, 6, 1), lty = c(1, 1, 3), bg = "white", lwd = 2)
dev.off()

# Some CIs:
est = average_gain(mq, spend)
round(est*100, 1)
round(100*(est[1] + c(-1.96, 1.96)*est[2]), 1)

est.civic = average_gain(mq.civic, spend)
round(est.civic*100, 1)
round(100*(est.civic[1] + c(-1.96, 1.96)*est.civic[2]), 1)

# paired BS, mq minus mq.self
est.diff = difference_gain(mq, mq.civic, spend)
round(100*est.diff, 1)
round(100*(est.diff[1] + c(-1.96, 1.96)*est.diff[2]), 1)

# MQ minus avg
est.avg.diff = difference_gain(mq, mq.avg, spend)
round(100*(est.avg.diff[1] + c(-1.96, 1.96)*est.avg.diff[2]), 1)

# Fraction treated in non-targeting and targeting
round(colMeans(predict(mq.avg, spend)), 2)
round(colMeans(predict(mq, spend)), 2)
