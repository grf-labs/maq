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
tau.true = predict(cf.train, data.truth$X, drop = TRUE)$predictions
cost.true = data.truth$cost
mq.true = maq(tau.true, cost.true, data.truth$tau)
gain.true = unlist(lapply(spends, function(s) average_gain(mq.true, s)[["estimate"]]))

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
    cost.eval = data.eval$cost
    mq = maq(tau.hat.eval, cost.eval, dr.eval, R = 200)

    est = unlist(lapply(spends, function(s) average_gain(mq, s)[["estimate"]]))
    se = unlist(lapply(spends, function(s) average_gain(mq, s)[["std.err"]]))

    coverage = as.integer(abs(est - gain.true) / se <= 1.96)
    bias = (est - gain.true)
    ci.length = se * 1.96 * 2

    res = c(res, list(data.frame(n = n, spend = spends, coverage, bias, ci.length)))
  }
}
Res = do.call(rbind, res)
res.df1 = aggregate(Res[, c(-1,-2)], by = list(n = Res$n, spend = Res$spend), FUN = mean)


# The mean coverage numbers, Table 1
xtable(xtabs(coverage ~ n + spend, res.df1))
