rm(list = ls())
set.seed(123)
source("generate_data.R")
library(maq)
library(grf)
library(xtable)

max.budget = 100
p = 10

# Generate "truth" test data
dgp = "map"
spends = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
data.truth = generate_data(100000, p, dgp = dgp)

# Fix a tau function
data.train = generate_data(10000, p, dgp = dgp)
cf.train = multi_arm_causal_forest(
  data.train$X,
  data.train$Y,
  data.train$W
)
# "population" quantities
cost.true = data.truth$cost
tau.true = predict(cf.train, data.truth$X, drop = TRUE)$predictions
mq.true = maq(tau.true, cost.true, max.budget, data.truth$tau)
mq1.true = maq(tau.true[, 1], cost.true[, 1], max.budget, data.truth$tau[, 1])
mq2.true = maq(tau.true[, 2], cost.true[, 2], max.budget, data.truth$tau[, 2])
mqr.true = maq(tau.true, cost.true, max.budget, data.truth$tau, target.with.covariates = FALSE)

gain.true = unlist(lapply(spends, function(s) average_gain(mq.true, s)[["estimate"]]))
gain1.true = unlist(lapply(spends, function(s) average_gain(mq1.true, s)[["estimate"]]))
gain2.true = unlist(lapply(spends, function(s) average_gain(mq2.true, s)[["estimate"]]))
gainr.true = unlist(lapply(spends, function(s) average_gain(mqr.true, s)[["estimate"]]))
true = c(gain.true - gainr.true, gain.true - gain1.true, gain.true - gain2.true, gain1.true - gainr.true, gain2.true - gainr.true)
name = rep(c("all_vs_r", "all_vs_1", "all_vs_2", "1_vs_r", "2_vs_r"), each = 10)

res = list()
for (n in c(1000, 2000, 5000, 10000)) {
  for (i in 1:1000) {
    data.eval = generate_data(n, p, dgp = dgp)
    tau.hat.eval = predict(cf.train, data.eval$X, drop = TRUE)$predictions

    cf.eval = multi_arm_causal_forest(
      data.eval$X,
      data.eval$Y,
      data.eval$W
    )
    dr.eval = get_scores(cf.eval, drop = TRUE)

    # fit maq
    cost.eval = data.eval$cost
    mq = maq(tau.hat.eval, cost.eval, max.budget, dr.eval, R = 200)
    mq1 = maq(tau.hat.eval[, 1], cost.eval[, 1], max.budget, dr.eval[, 1], R = 200)
    mq2 = maq(tau.hat.eval[, 2], cost.eval[, 2], max.budget, dr.eval[, 2], R = 200)
    mqr = maq(tau.hat.eval, cost.eval, max.budget, dr.eval, R = 200, target.with.covariates = FALSE)

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

xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "all_vs_r"))

xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "all_vs_1"))

xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "all_vs_2"))

xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "1_vs_r"))

xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "2_vs_r"))


# > xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "all_vs_r"))
# % latex table generated in R 3.5.1 by xtable 1.8-4 package
# % Mon May 29 02:53:43 2023
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrrr}
# \hline
# & 0.05 & 0.1 & 0.15 & 0.2 & 0.25 & 0.3 & 0.35 & 0.4 & 0.45 & 0.5 \\
# \hline
# 1000 & 0.95 & 0.96 & 0.96 & 0.97 & 0.96 & 0.95 & 0.96 & 0.96 & 0.96 & 0.96 \\
# 2000 & 0.95 & 0.95 & 0.95 & 0.94 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 \\
# 5000 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 \\
# 10000 & 0.95 & 0.94 & 0.94 & 0.94 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 \\
# \hline
# \end{tabular}
# \end{table}
# >
#   >
#   > xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "all_vs_1"))
# % latex table generated in R 3.5.1 by xtable 1.8-4 package
# % Mon May 29 02:53:43 2023
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrrr}
# \hline
# & 0.05 & 0.1 & 0.15 & 0.2 & 0.25 & 0.3 & 0.35 & 0.4 & 0.45 & 0.5 \\
# \hline
# 1000 & 0.95 & 0.95 & 0.96 & 0.96 & 0.95 & 0.95 & 0.94 & 0.95 & 0.95 & 0.95 \\
# 2000 & 0.95 & 0.95 & 0.95 & 0.93 & 0.95 & 0.95 & 0.94 & 0.94 & 0.95 & 0.95 \\
# 5000 & 0.95 & 0.96 & 0.95 & 0.95 & 0.94 & 0.95 & 0.95 & 0.96 & 0.96 & 0.95 \\
# 10000 & 0.95 & 0.94 & 0.94 & 0.94 & 0.95 & 0.95 & 0.96 & 0.96 & 0.95 & 0.95 \\
# \hline
# \end{tabular}
# \end{table}
# >
#   >
#   > xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "all_vs_2"))
# % latex table generated in R 3.5.1 by xtable 1.8-4 package
# % Mon May 29 02:53:43 2023
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrrr}
# \hline
# & 0.05 & 0.1 & 0.15 & 0.2 & 0.25 & 0.3 & 0.35 & 0.4 & 0.45 & 0.5 \\
# \hline
# 1000 & 0.95 & 0.95 & 0.96 & 0.95 & 0.96 & 0.95 & 0.95 & 0.96 & 0.94 & 0.94 \\
# 2000 & 0.96 & 0.95 & 0.97 & 0.96 & 0.96 & 0.96 & 0.96 & 0.95 & 0.95 & 0.96 \\
# 5000 & 0.96 & 0.95 & 0.95 & 0.95 & 0.95 & 0.94 & 0.95 & 0.95 & 0.96 & 0.95 \\
# 10000 & 0.95 & 0.95 & 0.96 & 0.96 & 0.95 & 0.95 & 0.96 & 0.96 & 0.95 & 0.95 \\
# \hline
# \end{tabular}
# \end{table}
# >
#   >
#   > xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "1_vs_r"))
# % latex table generated in R 3.5.1 by xtable 1.8-4 package
# % Mon May 29 02:53:43 2023
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrrr}
# \hline
# & 0.05 & 0.1 & 0.15 & 0.2 & 0.25 & 0.3 & 0.35 & 0.4 & 0.45 & 0.5 \\
# \hline
# 1000 & 0.96 & 0.95 & 0.95 & 0.96 & 0.95 & 0.96 & 0.96 & 0.96 & 0.96 & 0.95 \\
# 2000 & 0.95 & 0.96 & 0.95 & 0.96 & 0.96 & 0.96 & 0.96 & 0.96 & 0.96 & 0.95 \\
# 5000 & 0.95 & 0.95 & 0.95 & 0.94 & 0.94 & 0.95 & 0.95 & 0.95 & 0.96 & 0.95 \\
# 10000 & 0.95 & 0.95 & 0.95 & 0.94 & 0.94 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 \\
# \hline
# \end{tabular}
# \end{table}
# >
#   >
#   > xtable(xtabs(coverage ~ n + spend, res.df2, subset = res.df2$name == "2_vs_r"))
# % latex table generated in R 3.5.1 by xtable 1.8-4 package
# % Mon May 29 02:53:43 2023
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrrrrrrrr}
# \hline
# & 0.05 & 0.1 & 0.15 & 0.2 & 0.25 & 0.3 & 0.35 & 0.4 & 0.45 & 0.5 \\
# \hline
# 1000 & 0.95 & 0.95 & 0.96 & 0.95 & 0.95 & 0.95 & 0.96 & 0.96 & 0.95 & 0.96 \\
# 2000 & 0.95 & 0.96 & 0.95 & 0.94 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 & 0.95 \\
# 5000 & 0.95 & 0.95 & 0.95 & 0.95 & 0.94 & 0.95 & 0.95 & 0.95 & 0.94 & 0.95 \\
# 10000 & 0.94 & 0.95 & 0.94 & 0.95 & 0.96 & 0.95 & 0.96 & 0.94 & 0.95 & 0.95 \\
# \hline
# \end{tabular}
# \end{table}



