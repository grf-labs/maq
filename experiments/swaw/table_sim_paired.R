rm(list = ls())
set.seed(123)

library(maq)
library(grf)
library(xtable)
source("generate_data.R")

# Number of Monte Carlo repetitions
num.mc = 1000 # This takes a few days to complete

# Generate "ground truth" data
spends = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
data.truth = generate_data(n = 100000, p = 10, dgp = "map")

# Fix a tau function
data.train = generate_data(n = 10000, p = 10, dgp = "map")
cf.train = multi_arm_causal_forest(
  data.train$X,
  data.train$Y,
  data.train$W
)
# "population" quantities
cost.true = data.truth$cost
tau.true = predict(cf.train, data.truth$X, drop = TRUE)$predictions
mq.true = maq(tau.true, cost.true, data.truth$tau)
mq1.true = maq(tau.true[, 1], cost.true[, 1], data.truth$tau[, 1])
mq2.true = maq(tau.true[, 2], cost.true[, 2], data.truth$tau[, 2])
mqr.true = maq(tau.true, cost.true, data.truth$tau, target.with.covariates = FALSE)

gain.true = unlist(lapply(spends, function(s) average_gain(mq.true, s)[["estimate"]]))
gain1.true = unlist(lapply(spends, function(s) average_gain(mq1.true, s)[["estimate"]]))
gain2.true = unlist(lapply(spends, function(s) average_gain(mq2.true, s)[["estimate"]]))
gainr.true = unlist(lapply(spends, function(s) average_gain(mqr.true, s)[["estimate"]]))
true = c(gain.true - gainr.true, gain.true - gain1.true, gain.true - gain2.true, gain1.true - gainr.true, gain2.true - gainr.true)
name = rep(c("all_vs_r", "all_vs_1", "all_vs_2", "1_vs_r", "2_vs_r"), each = 10)

res = list()
for (n in c(1000, 2000, 5000, 10000)) {
  for (i in 1:num.mc) {
    data.eval = generate_data(n, p = 10, dgp = "map")
    tau.hat.eval = predict(cf.train, data.eval$X, drop = TRUE)$predictions

    cf.eval = multi_arm_causal_forest(
      data.eval$X,
      data.eval$Y,
      data.eval$W
    )
    dr.eval = get_scores(cf.eval, drop = TRUE)

    # fit maq
    cost.eval = data.eval$cost
    mq = maq(tau.hat.eval, cost.eval, dr.eval, R = 200)
    mq1 = maq(tau.hat.eval[, 1], cost.eval[, 1], dr.eval[, 1], R = 200)
    mq2 = maq(tau.hat.eval[, 2], cost.eval[, 2], dr.eval[, 2], R = 200)
    mqr = maq(tau.hat.eval, cost.eval, dr.eval, R = 200, target.with.covariates = FALSE)

    est.all.r = unlist(lapply(spends, function(s) difference_gain(mq, mqr, s)[["estimate"]]))
    se.all.r = unlist(lapply(spends, function(s) difference_gain(mq, mqr, s)[["std.err"]]))

    est.all.1.r = unlist(lapply(spends, function(s) difference_gain(mq, mq1, s)[["estimate"]]))
    se.all.1.r = unlist(lapply(spends, function(s) difference_gain(mq, mq1, s)[["std.err"]]))

    est.all.2.r = unlist(lapply(spends, function(s) difference_gain(mq, mq2, s)[["estimate"]]))
    se.all.2.r = unlist(lapply(spends, function(s) difference_gain(mq, mq2, s)[["std.err"]]))

    est.1.r = unlist(lapply(spends, function(s) difference_gain(mq1, mqr, s)[["estimate"]]))
    se.1.r = unlist(lapply(spends, function(s) difference_gain(mq1, mqr, s)[["std.err"]]))

    est.2.r = unlist(lapply(spends, function(s) difference_gain(mq2, mqr, s)[["estimate"]]))
    se.2.r = unlist(lapply(spends, function(s) difference_gain(mq2, mqr, s)[["std.err"]]))

    est = c(est.all.r, est.all.1.r, est.all.2.r, est.1.r, est.2.r)
    se = c(se.all.r, se.all.1.r, se.all.2.r, se.1.r, se.2.r)

    coverage = as.integer(abs(est - true) / se <= 1.96)
    bias = (est - true)
    ci.length = se * 1.96 * 2

    res = c(res, list(data.frame(n = n, spend = spends, name = name, coverage, bias, ci.length)))
  }
}
Res = do.call(rbind, res)
res.df2 = aggregate(Res[, c(-1,-2, -3)], by = list(n = Res$n, spend = Res$spend, name = Res$name), FUN = mean)

# The mean coverage numbers, Table 3

# Panel A
xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "all_vs_r"))

# Panel B
xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "all_vs_1"))

# Panel C
xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "all_vs_2"))

# Panel D
xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "1_vs_r"))

# Panel E
xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "2_vs_r"))
