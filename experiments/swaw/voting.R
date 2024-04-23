rm(list = ls())
set.seed(123)
library(maq)
library(grf)

# Voting data file "D001F01" from https://isps.yale.edu/research/data/d001
# Alan S Gerber, Donald P Green, and Christopher W Larimer. Social pressure and voter turnout:
# Evidence from a large-scale field experiment. American Political Science Review, 102(1):33–48, 2008.
raw = read.csv("GerberGreenLarimer_APSR_2008_social_pressure.csv.gz")
raw$treatment = trimws(raw$treatment)

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

# Draw a random sample of households for train/evaluation (50/50)
samples.by.hh = split(seq_along(cluster), cluster)
num.hh = length(samples.by.hh)
train = unlist(samples.by.hh[sample(1:num.hh, num.hh / 2)], use.names = FALSE)
eval = -train

# Specify the propensity scores (the trial's randomization probabilities)
W.hat = c(5/9, 1/9, 1/9, 1/9, 1/9)
# Train a CATE function
cf.train = multi_arm_causal_forest(X[train, ],
                                   Y[train],
                                   W[train],
                                   W.hat = W.hat,
                                   clusters = cluster[train],
                                   num.trees = 500)
# Predict CATEs on evaluation set
tau.hat.eval = predict(cf.train, X[eval, ], drop = TRUE)$predictions

# Get Qini evaluation scores using inverse-propensity weighting (IPW)
Y.k.ipw.eval = get_ipw_scores(Y[eval], W[eval], W.hat)

# Specify the cost of each arm.
cost = c(1, 15, 30, 45)

# Fit a Qini curve using all arms
mq = maq(tau.hat.eval, cost, Y.k.ipw.eval, clusters = cluster[eval], R = 200)

# Fit a Qini curve using only arm 1 ("Civic")
mq.civic = maq(tau.hat.eval[, 1], cost[1], Y.k.ipw.eval[, 1], clusters = cluster[eval], R = 200)

# Fit a Qini curve using only arm 2 ("Hawthorne")
mq.hawthorne = maq(tau.hat.eval[, 2], cost[2], Y.k.ipw.eval[, 2], clusters = cluster[eval], R = 200)

# Fit a Qini curve using only arm 3 ("Self")
mq.self = maq(tau.hat.eval[, 3], cost[3], Y.k.ipw.eval[, 3], clusters = cluster[eval], R = 200)

# Fit a Qini curve using only arm 4 ("Neighbors")
mq.neighbors = maq(tau.hat.eval[, 4], cost[4], Y.k.ipw.eval[, 4], clusters = cluster[eval], R = 200)

# Fit a non-targeting baseline
mq.avg = maq(tau.hat.eval, cost, Y.k.ipw.eval, clusters = cluster[eval], R = 200, target.with.covariates = FALSE)


spend = 5
# Plot all the estimated curves
pdf("voting_figure_a.pdf", pointsize = 18)
plot(mq, xlab = "Intrusion cost", ylab = "Increase in probability of voting",
     ci.args = NULL,  ylim = c(0, 0.06), lwd = 3, lty = 1)
plot(mq.civic, add = TRUE, col = 2, ci.args = NULL, lwd = 3, lty = 2)
plot(mq.hawthorne, add = TRUE, col = 3, ci.args = NULL, lwd = 3, lty = 3)
plot(mq.self, add = TRUE, col = 4, ci.args = NULL, lwd = 3, lty = 4)
plot(mq.neighbors, add = TRUE, col = 7, ci.args = NULL, lwd = 3, lty = 5)
abline(v = spend, lty = 2)
legend("topleft",
       c("All", "Civic", "Hawthorne", "Self", "Neighbors"),
       col = c(1,2,3,4,7), lty = 1:5, bg = "white", lwd = 2)
dev.off()

# Plot the multi-armed curves
pdf("voting_figure_b.pdf", pointsize = 18)
plot(mq, xlab = "Intrusion cost", ylab = "Increase in probability of voting",
     ylim = c(0, 0.06), lwd = 3)
plot(mq.avg, add = TRUE, col = 6, lwd = 3, lty = 2)
abline(v = spend, lty = 2)
legend("topleft",
       c("All (targeting)", "All (no targeting)"),
       col = c(1, 6, 1), lty = c(1, 2, 3), lwd = 2, bg = "white")
dev.off()

# Some CIs:
# Multi-armed targeting at B = 5
est = average_gain(mq, spend)
round(est*100, 1)
round(100*(est[1] + c(-1.96, 1.96)*est[2]), 1)

# Targeting only with "Civic" at B = 5
est.civic = average_gain(mq.civic, spend)
round(est.civic*100, 1)
round(100*(est.civic[1] + c(-1.96, 1.96)*est.civic[2]), 1)

# Get a paired difference "All - Civic" at B = 5
est.diff = difference_gain(mq, mq.civic, spend)
round(100*est.diff, 1)
round(100*(est.diff[1] + c(-1.96, 1.96)*est.diff[2]), 1)

# Multi-armed targeting vs baseline at B = 5
est.avg.diff = difference_gain(mq, mq.avg, spend)
round(100*(est.avg.diff[1] + c(-1.96, 1.96)*est.avg.diff[2]), 1)

# Fraction treated in non-targeting and targeting
round(colMeans(predict(mq.avg, spend)), 2)
round(colMeans(predict(mq, spend)), 2)
