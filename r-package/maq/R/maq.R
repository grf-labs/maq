#' Fit a Multi-Action Qini.
#'
#'
#' @param reward A matrix of reward estimates.
#' @param cost A matrix of cost estimates. If the costs are the same for each unit, then this can also
#'  be a ncol(reward)-length vector.
#' @param budget The maximum spend per unit to fit the MAQ path on.
#' @param reward.scores A matrix of rewards to evaluate the MAQ on. For valid statistical inference, the
#'  reward and cost estimates should be obtained independently from this evaluation data.
#' @param target.with.covariates If TRUE, then the optimal policy takes covariates into
#'  account. If FALSE, then the optimal policy only takes the average reward and cost into account when
#'  allocating treatment. Default is TRUE.
#' @param R Number of bootstrap replicates for computing standard errors. Default is 0
#'  (only point estimates are computed).
#' @param paired.inference Whether to allow for paired tests with other cost curves. If TRUE (Default)
#'  then the path of bootstrap replicates are stored in order to perform paired comparisons.
#'  This takes memory on the order of O(RnK) and requires the comparison objects to be fit with the
#'  same seed and R values as well as the same number of samples.
#' @param sample.weights Weights given to an observation in estimation.
#'  If NULL, each observation is given the same weight. Default is NULL.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#'  Default is NULL (ignored).
#' @param tie.breaker An optional permutation of the the integers 1 to nrow(rewards) used to
#'  break potential ties in the optimal treatment allocation. If NULL, the ties are broken by
#'  the lowest sample id (i.e. the sample appearing first in the data). Default is NULL.
#' @param num.threads Number of threads used in bootstrap replicates. By default, the number of threads
#'  is set to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator. Default is 42.
#'
#' @return A fit maq object.
#'
#' @examples
#' \donttest{
#' if (require("grf", quietly = TRUE)) {
#'
#' # Fit a CATE estimator (using GRF) on a training sample.
#' n <- 3000
#' p <- 5
#' X <- matrix(runif(n * p), n, p)
#' W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
#' Y <- X[, 1] + X[, 2] * (W == "B") + 1.5 * X[, 3] * (W == "C") + rnorm(n)
#' train <- sample(1:n, n/2)
#' eval <- -train
#'
#' tau.forest <- grf::multi_arm_causal_forest(X[train, ], Y[train], W[train])
#'
#' # Predict CATEs on held out evaluation data.
#' tau.hat <- predict(tau.forest, X[eval, ], drop = TRUE)$predictions
#'
#' # Form cost estimates - the following are a toy example.
#' cost.hat <- cbind(X[eval, 4] / 4, X[eval, 5])
#'
#' # Fit an evaluation forest to compute doubly robust evaluation set scores.
#' eval.forest <- grf::multi_arm_causal_forest(X[eval, ], Y[eval], W[eval])
#' DR.scores <- grf::get_scores(eval.forest, drop = TRUE)
#'
#' # Fit a MAQ with evaluation set estimates using 200 bootstrap replicates for confidence intervals.
#' max.budget <- 1
#' mq <- maq(tau.hat, cost.hat, max.budget, DR.scores, R = 200)
#'
#' # Plot the MAQ curve.
#' plot(mq)
#' legend("topleft", c("All arms", "95% CI"), lty = c(1, 3))
#'
#' # Get an estimate of optimal reward at a given spend per unit along with standard errors.
#' average_gain(mq, spend = 0.2)
#'
#' # Get the optimal treatment allocation matrix at a given spend per unit.
#' pi.mat <- predict(mq, spend = 0.2)
#'
#' # If the treatment randomization probabilities are known, then an alternative to
#' # evaluation via AIPW scores is to use inverse-propensity weighting (IPW).
#' W.hat.true <- rep(1/3, 3)
#' observed.W <- match(W, levels(W))
#' Y.k.mat <- matrix(0, length(W), nlevels(W))
#' Y.k.mat[cbind(seq_along(observed.W), observed.W)] <- Y
#' Y.k.ipw <- sweep(Y.k.mat, 2, W.hat.true, "/")
#' Y.k.ipw.eval <- Y.k.ipw[eval, -1] - Y.k.ipw[eval, 1]
#'
#' mq.ipw <- maq(tau.hat, cost.hat, max.budget, Y.k.ipw.eval)
#' plot(mq.ipw, add = TRUE, col = 2)
#' legend("topleft", c("All arms", "95% CI", "All arms (IPW)"), col = c(1, 1, 2), lty = c(1, 3, 1))
#'
#' # Estimate some baseline policies.
#' # a) A policy that ignores covariates and only only takes the average reward/cost into account.
#' mq.avg <- maq(tau.hat, cost.hat, max.budget, DR.scores, target.with.covariates = FALSE, R = 200)
#'
#' # b) A policy that only use arm 1.
#' mq.arm1 <- maq(tau.hat[, 1], cost.hat[, 1], max.budget, DR.scores[, 1], R = 200)
#'
#' # c) A policy that only use arm 2.
#' mq.arm2 <- maq(tau.hat[, 2], cost.hat[, 2], max.budget, DR.scores[, 2], R = 200)
#'
#' plot(mq, ci.args = NULL)
#' plot(mq.avg, col = 2, add = TRUE, ci.args = NULL)
#' plot(mq.arm1, col = 3, add = TRUE, ci.args = NULL)
#' plot(mq.arm2, col = 4, add = TRUE, ci.args = NULL)
#' legend("topleft", c("All arms (targeting)", "All arms (without targeting)", "Arm 1", "Arm 2"),
#'        col = 1:4, lty = 1)
#'
#' # Estimate the value of employing all arms over a random allocation.
#' difference_gain(mq, mq.avg, spend = 0.2)
#'
#' # Estimate the value of adding arm 1 to the optimal policy mix.
#' difference_gain(mq, mq.arm1, spend = 0.2)
#'
#' # Estimate the value of adding arm 2 to the optimal policy mix.
#' difference_gain(mq, mq.arm2, spend = 0.2)
#'
#' }
#' }
#'
#' @export
maq <- function(reward,
                cost,
                budget,
                reward.scores,
                target.with.covariates = TRUE,
                R = 0,
                paired.inference = TRUE,
                sample.weights = NULL,
                clusters = NULL,
                tie.breaker = NULL,
                num.threads = NULL,
                seed = 42) {
  if (is.vector(cost) && length(cost) == NCOL(reward.scores)) {
    cost <- matrix(cost, NROW(reward.scores), NCOL(reward.scores), byrow = TRUE)
  }

  if (NROW(reward) != NROW(cost) || NCOL(reward) != NCOL(cost)
        || NROW(reward) != NROW(reward.scores) || NCOL(reward) != NCOL(reward.scores)
        || anyNA(reward) || anyNA(cost) || anyNA(reward.scores)) {
    stop("rewards and costs should be matrices of equal size with no missing values.")
  }

  if (any(cost <= 0)) {
    stop("Costs should be > 0.")
  }

  if (R < 0) {
    stop("The number of bootstrap replicates R should be a non-negative integer.")
  }

  if (is.null(sample.weights)) {
    sample.weights <- vector(mode = "numeric", length = 0)
  } else if (length(sample.weights) != NROW(reward) || anyNA(sample.weights)
               || any(sample.weights <= 0)) {
    stop("sample.weights should have length=nrow(reward) and be non-missing and positive.")
  } else {
    sample.weights <- sample.weights / sum(sample.weights)
  }

  if (is.null(clusters)) {
    clusters <- vector(mode = "numeric", length = 0)
  } else {
    if (mode(clusters) != "numeric") {
      stop("clusters must be able to be coerced to a numeric vector.")
    }
    clusters <- as.numeric(clusters)
    if (!all(clusters == floor(clusters))) {
      stop("clusters vector cannot contain floating point values.")
    } else if (length(clusters) != NROW(reward)) {
      stop("clusters vector has incorrect length.")
    } else {
      # convert to integers between 0 and n clusters
      clusters <- as.numeric(as.factor(clusters)) - 1
    }
    cluster.size.counts <- table(clusters)
    if (floor(length(cluster.size.counts) / 2) <= 1) {
      stop("Cannot bootstrap sample with only one effective unit.")
    }
  }

  if (is.null(tie.breaker)) {
    tie.breaker <- vector(mode = "integer", length = 0)
  } else if (length(tie.breaker) != NROW(reward)) {
    stop("tie.breaker should have length=nrow(reward).")
  }

  if (is.null(num.threads)) {
    num.threads <- 0
  } else if (num.threads < 0) {
    stop("num.threads should be a non-negative integer.")
  }

  if (!is.numeric(seed) || seed < 0) {
    stop("seed should be a non-negative integer.")
  }

  ret <- solver_rcpp(as.matrix(reward), as.matrix(reward.scores), as.matrix(cost),
                     sample.weights, tie.breaker, clusters,
                     budget, target.with.covariates, paired.inference, R, num.threads, seed)

  output <- list()
  class(output) <- "maq"
  output[["_path"]] <- ret
  output[["seed"]] <- seed
  output[["target.with.covariates"]] <- target.with.covariates
  output[["paired.inference"]] <- paired.inference
  output[["R"]] <- R
  output[["dim"]] <- c(NROW(cost), NCOL(cost))
  output[["budget"]] <- budget

  output
}

#' Predict optimal treatment allocation.
#'
#'
#' @param object A maq object.
#' @param spend The spend level.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A sparse matrix.
#' @method predict maq
#' @export
predict.maq <- function(object,
                        spend,
                        ...) {
  if (!object[["_path"]]$complete.path && spend > object$budget) {
    stop("maq path is not fit beyond given spend level.")
  }
  spend.grid <- object[["_path"]]$spend
  path.idx <- findInterval(spend, spend.grid) # nearest path index (lower bound)
  if (path.idx == 0) {
    return (Matrix::sparseMatrix(i = NULL, j = NULL, x = 0, dims = object[["dim"]]))
  }

  ipath <- object[["_path"]]$ipath[1:path.idx] + 1 # +1: R index.
  kpath <- object[["_path"]]$kpath[1:path.idx] + 1
  ix <- !duplicated(ipath, fromLast = TRUE)
  pi.mat <- Matrix::sparseMatrix(ipath[ix], kpath[ix], x = 1, dims = object[["dim"]])
  if (path.idx == length(spend.grid)) {
    return (pi.mat)
  }
  # fractional adjustment?
  spend.diff <- spend - spend.grid[path.idx]
  next.unit <- object[["_path"]]$ipath[path.idx + 1] + 1
  next.arm <- object[["_path"]]$kpath[path.idx + 1] + 1
  prev.arm <- Matrix::which(pi.mat[next.unit, ] == 1) # already assigned?

  fraction <- spend.diff / (spend.grid[path.idx + 1] - spend.grid[path.idx])
  pi.mat[next.unit, next.arm] <- fraction
  if (length(prev.arm) > 0) {
    pi.mat[next.unit, prev.arm] <- 1 - fraction
  }

  pi.mat
}

#' Get estimate of gain given a spend level.
#'
#'
#' @param object A maq object.
#' @param spend The spend level.
#'
#' @return An estimate of average gain along with standard errors.
#' @export
average_gain <- function(object,
                         spend) {
  if (!object[["_path"]]$complete.path && spend > object$budget) {
    stop("maq path is not fit beyond given spend level.")
  }
  spend.grid <- object[["_path"]]$spend
  path.idx <- findInterval(spend, spend.grid) # nearest path index (lower bound)

  gain.path <- object[["_path"]]$gain
  se.path <- object[["_path"]]$std.err
  if (path.idx == 0) {
    estimate <- 0
    std.err <- 0
  } else if (path.idx == length(spend.grid)) {
    estimate <- gain.path[path.idx]
    std.err <- se.path[path.idx]
  } else {
    interp.ratio <- (spend - spend.grid[path.idx]) / (spend.grid[path.idx + 1] - spend.grid[path.idx])
    estimate <- gain.path[path.idx] + (gain.path[path.idx + 1] - gain.path[path.idx]) * interp.ratio
    std.err <- se.path[path.idx] + (se.path[path.idx + 1] - se.path[path.idx]) * interp.ratio
  }

  c(estimate = estimate, std.err = std.err)
}

#' Get estimate of difference in gain given a spend level.
#'
#'
#' @param object.lhs A maq object to subtract from.
#' @param object.rhs A maq object to subtract with.
#' @param spend The spend level.
#'
#' @return An estimate of difference in gain along with standard errors.
#' @export
difference_gain <- function(object.lhs,
                            object.rhs,
                            spend) {
  if (!object.lhs[["_path"]]$complete.path && spend > object.lhs$budget) {
    stop("lhs maq path is not fit beyond given spend level.")
  }
  if (!object.rhs[["_path"]]$complete.path && spend > object.rhs$budget) {
    stop("rhs maq path is not fit beyond given spend level.")
  }
  if (object.lhs[["seed"]] != object.rhs[["seed"]] ||
      object.lhs[["R"]] != object.rhs[["R"]] ||
      object.lhs[["dim"]][[1]] != object.rhs[["dim"]][[1]] ||
      !object.lhs[["paired.inference"]] ||
      !object.rhs[["paired.inference"]]) {
    stop(paste("Paired comparisons require maq objects to be fit with paired.inference=TRUE",
               "as well as with the same random seed, bootstrap replicates, and data"))
  }

  estimate <- average_gain(object.lhs, spend)[[1]] - average_gain(object.rhs, spend)[[1]]
  # Compute paired std.errors
  .get_estimates <- function(object) {
    gain.bs <- object[["_path"]]$gain.bs
    spend.grid <- object[["_path"]]$spend
    path.idx <- findInterval(spend, spend.grid) # nearest path index (lower bound)
    if (path.idx == 0) {
      estimates <- 0
    } else if (path.idx == length(spend.grid)) {
      estimates <- unlist(lapply(gain.bs, function(gain.path.bs) gain.path.bs[path.idx]))
    } else {
      interp.ratio <- (spend - spend.grid[path.idx]) / (spend.grid[path.idx + 1] - spend.grid[path.idx])
      estimates <- unlist(lapply(gain.bs, function(gain.path.bs) {
        gain.path.bs[path.idx] + (gain.path.bs[path.idx + 1] - gain.path.bs[path.idx]) * interp.ratio
      }))
    }

    estimates
  }
  estimates.lhs <- .get_estimates(object.lhs)
  estimates.rhs <- .get_estimates(object.rhs)
  std.err <- stats::sd(estimates.lhs - estimates.rhs, na.rm = TRUE)
  if (is.na(std.err)) {
    std.err <- 0
  }

  c(estimate = estimate, std.err = std.err)
}

#' MAQ Summary.
#' @param object A maq object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A data.frame with the computed path.
#' @method summary maq
#' @export
summary.maq <- function(object,
                        ...) {

  data.frame(
    spend = object[["_path"]]$spend,
    gain = object[["_path"]]$gain,
    std.err = object[["_path"]]$std.err,
    unit.allocation = object[["_path"]]$ipath + 1,
    arm.allocation = object[["_path"]]$kpath + 1
  )
}

#' Print a maq object.
#' @param x A maq object.
#' @param ... Additional arguments (currently ignored).
#'
#' @method print maq
#' @export
print.maq <- function(x,
                      ...) {

  cat("MAQ object fit on", x$dim[1], "units and", x$dim[2], "arms with max budget", x$budget)
}

#' Plot the gain/spend curve.
#' @param x A maq object.
#' @param ... Additional arguments passed to plot.
#' @param add Whether to add to an already existing plot. Default is FALSE.
#' @param horizontal.line Whether to draw a horizontal line where the cost curve continues.
#'  Only applies if add = TRUE. Default is TRUE.
#' @param ci.args A list of optional arguments to lines() for drawing 95 % confidence bars.
#'  Set to NULL to ignore CIs.
#' @param grid.step The grid increment size to plot the curve on. Default is
#'  max(floor(length(path.length) / 1000), 1).
#'
#' @method plot maq
#' @export
plot.maq <- function(x,
                     ...,
                     add = FALSE,
                     horizontal.line = TRUE,
                     ci.args = list(),
                     grid.step = NULL
                     ) {
  spend.grid <- x[["_path"]]$spend
  gain.grid <- x[["_path"]]$gain
  std.err.grid <- x[["_path"]]$std.err
  if (length(spend.grid) < 1) {
    return(invisible(x))
  }

  if (is.null(grid.step)) {
    grid.step <- max(floor(length(spend.grid) / 1000), 1)
  }
  plot.grid <- seq(1, length(spend.grid), by = grid.step)
  spend <- spend.grid[plot.grid]
  gain <- gain.grid[plot.grid]
  std.err <- std.err.grid[plot.grid]
  if (add && horizontal.line) {
    len <- length(spend)
    xmax <- graphics::par("usr")[2]
    spend <- c(spend, seq(spend[len], xmax, length.out = 100))
    gain <- c(gain, rep(gain[len], 100))
    std.err <- c(std.err, rep(std.err[len], 100))
  }
  lb <- gain - 1.96 * std.err
  ub <- gain + 1.96 * std.err

  plot.args <- list(type = "l", ylim = c(min(lb), max(ub)), xlab = "spend", ylab = "gain", col = 1)
  new.args <- list(...)
  plot.args[names(new.args)] <- new.args

  lines.args <- list(lty = 3, col = plot.args$col)
  lines.args[names(ci.args)] <- ci.args

  if (!add || grDevices::dev.cur() == 1L) {
    do.call(plot, c(list(x = spend, y = gain), plot.args))
  } else {
    do.call(graphics::lines, c(list(x = spend, y = gain), plot.args))
  }

  if (!is.null(ci.args)) {
    do.call(graphics::lines, c(list(x = spend, y = lb), lines.args))
    do.call(graphics::lines, c(list(x = spend, y = ub), lines.args))
  }
}
