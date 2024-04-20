rm(list = ls())
set.seed(123)

library(maq)
library(grf)
library(xtable)
source("generate_data.R")

# Number of Monte Carlo repetitions
num.mc = 1000 # This takes a few days to complete

# XGBoost with cv-tuned hyperparameters.
# (adapted from Xinkun Nie's R-learner https://github.com/xnie/rlearner/blob/master/R/cvboost.R)
cvboost.fit.predict = function(
    X,
    Y,
    W,
    k_folds = 3,
    objective ="reg:squarederror",
    ntrees_max = 500,
    num_search_rounds = 10,
    print_every_n = 100,
    early_stopping_rounds = 10) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Requires the package 'xgboost'.")
  }
  eval = "rmse"
  if (objective == "reg:squarederror") {
    dtrain <- xgboost::xgb.DMatrix(data = X, label = Y)
  } else {
    dtrain <- xgboost::xgb.DMatrix(data = X, label = as.integer(Y == 0))
  }
  best_param = list()
  best_loss = Inf
  folds = split(1:length(Y), sample(k_folds, length(Y), replace = TRUE))

  for (iter in 1:num_search_rounds) {
    param <- list(objective = objective,
                  eval_metric = eval,
                  subsample = sample(c(0.5, 0.75, 1), 1),
                  colsample_bytree = sample(c(0.6, 0.8, 1), 1),
                  eta = sample(c(5e-3, 1e-2, 0.015, 0.025, 5e-2, 8e-2, 1e-1, 2e-1), 1),
                  max_depth = sample(c(3:20), 1),
                  gamma = runif(1, 0.0, 0.2),
                  min_child_weight = sample(1:20, 1),
                  max_delta_step = sample(1:10, 1))

    xgb_cv_args = list(data = dtrain,
                       param = param,
                       missing = NA,
                       folds = folds,
                       prediction = TRUE,
                       early_stopping_rounds = early_stopping_rounds,
                       maximize = FALSE,
                       nrounds = ntrees_max,
                       print_every_n = print_every_n,
                       nthread = parallel::detectCores(),
                       callbacks = list(xgboost::cb.cv.predict(save_models = TRUE)))

    xgb_cvfit <- do.call(xgboost::xgb.cv, xgb_cv_args)

    metric = paste('test_', eval, '_mean', sep='')
    min_loss = min(xgb_cvfit$evaluation_log[, ..metric])

    if (min_loss < best_loss) {
      best_loss = min_loss
      best_param = param
      best_xgb_cvfit = xgb_cvfit
    }
  }
  xgb_train_args = list(data = dtrain,
                        params = best_param,
                        nrounds = best_xgb_cvfit$best_ntreelimit)

  # Fit xgboost models for each arm.
  if (objective == "reg:squarederror") {
    preds = lapply(0:2, function(k) {
      idx = W == k
      data = xgboost::xgb.DMatrix(data = X[idx, ], label = Y[idx])
      xgb_train_args$data = data
      fit = do.call(xgboost::xgb.train, xgb_train_args)
      predict(fit, X)
    })
  } else {
    preds = lapply(0:2, function(k) {
      idx = W == k
      data = xgboost::xgb.DMatrix(data = X[idx, ], label = as.integer(Y[idx] == k))
      xgb_train_args$data = data
      fit = do.call(xgboost::xgb.train, xgb_train_args)
      predict(fit, X)
    })
  }

  out = matrix(unlist(preds), nrow(X), 3)
  if (objective == "binary:logistic") {
    out = sweep(out, 1, rowSums(out), "/")
  }

  out
}

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

    muhat = cvboost.fit.predict(data.eval$X, data.eval$Y, data.eval$W)
    What = cvboost.fit.predict(data.eval$X, data.eval$W, data.eval$W, objective = "binary:logistic")

    dr.eval = get_aipw_scores(data.eval$Y, data.eval$W, muhat, What)
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

# Mean coverage, Table C1
xtable(xtabs(coverage ~ n + spend, res.df1))
