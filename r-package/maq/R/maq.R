#' Fit a Multi-Action QINI.
#'
#'
#' @param reward A matrix of reward estimates.
#' @param cost A matrix of cost estimates.
#' @param budget The maximum spend/unit to fit the MAQ path on.
#' @param reward.scores A matrix of evaluation set reward score estimates.
#' @param R Number of bootstrap replicates for SEs. Default is 200.
#' @param sample.weights Weights given to an observation in estimation.
#'  If NULL, each observation is given the same weight. Default is NULL.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to.
#'  Default is NULL (ignored).
#' @param tie.breaker An optional permutation of the the integers 1 to nrow(rewards) used to
#'  break potential ties in the optimal treatment allocation. If NULL, the ties are broken by
#'  the lowest sample id (i.e. the sample appearing first in the data). Default is NULL.
#' @param num.threads Number of threads used in bootstrap replicates. By default, the number of threads
#'  is set to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator.
#'
#' @return A fit maq object.
#'
#' @examples
#' \donttest{
#' if (require("grf", quietly = TRUE)) {
#' # Fit a CATE estimator (using GRF) on a training sample.
#' n <- 2000
#' p <- 5
#' X <- matrix(runif(n * p), n, p)
#' W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
#' Y <- X[, 1] + X[, 2] * (W == "B") + X[, 3] * (W == "C") + rnorm(n)
#' train <- sample(1:n, n/2)
#' eval <- -train
#'
#' tau.forest <- grf::multi_arm_causal_forest(X[train, ], Y[train], W[train])
#'
#' # Predict CATEs on held out evaluation data.
#' tau.hat <- predict(tau.forest, X[eval, ])$predictions[,,]
#'
#' # Form cost estimates - the following are a toy example.
#' cost.hat <- X[eval, 4:5]
#'
#' # Fit a evaluation forest to compute doubly robust evaluation set scores.
#' eval.forest <- grf::multi_arm_causal_forest(X[eval, ], Y[eval], W[eval])
#' DR.scores <- grf::get_scores(eval.forest)[,,]
#'
#' # Fit a MAQ using evaluation set estimates.
#' max.budget <- 1
#' mq <- maq(tau.hat, cost.hat, max.budget, DR.scores)
#'
#' # Plot the MAQ curve.
#' plot(mq)
#'
#' # Get an estimate of optimal reward at a given spend per unit along with standard errors.
#' average_gain(mq, spend = 0.3)
#'
#' # Get the optimal treatment allocation matrix at a given spend per unit.
#' pi.mat <- predict(mq, spend = 0.3)
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
#' }
#' }
#'
#' @export
maq <- function(reward,
                cost,
                budget,
                reward.scores,
                R = 200,
                sample.weights = NULL,
                clusters = NULL,
                tie.breaker = NULL,
                num.threads = NULL,
                seed = runif(1, 0, .Machine$integer.max)) {
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
    samples.per.cluster <- 0
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
    samples.per.cluster <- max(cluster.size.counts)
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
                     samples.per.cluster, budget, R, num.threads, seed)

  output <- list()
  class(output) <- "maq"
  output[["_path"]] <- ret
  output[["seed"]] <- seed
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

#' MAQ Summary.
#' @param object The maq object.
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
    allocated.unit = object[["_path"]]$ipath + 1,
    allocated.arm = object[["_path"]]$kpath + 1
  )
}

#' Print a maq object.
#' @param x The maq object.
#' @param ... Additional arguments (currently ignored).
#'
#' @method print maq
#' @export
print.maq <- function(x,
                      ...) {

  cat("MAQ object fit on", x$dim[1], "units and", x$dim[2], "arms with max budget", x$budget)
}

#' Plot the gain/spend curve.
#' @param x The output of maq.
#' @param ... Additional arguments passed to plot.
#'
#' @method plot maq
#' @export
plot.maq <- function(x,
                     ...) {
  spend.grid <- x[["_path"]]$spend
  gain.grid <- x[["_path"]]$gain
  std.err.grid <- x[["_path"]]$std.err
  if (length(spend.grid) < 1) {
    return(invisible(x))
  }

  plot.grid <- seq(1, length(spend.grid), by = max(floor(length(spend.grid) / 1000), 1))
  spend <- spend.grid[plot.grid]
  gain <- gain.grid[plot.grid]
  std.err <- std.err.grid[plot.grid]
  lb <- gain - 1.96 * std.err
  ub <- gain + 1.96 * std.err

  plot(spend, gain, type = "l", ylim = c(min(lb), max(ub)), ...)
  graphics::lines(spend, lb, lty = 3)
  graphics::lines(spend, ub, lty = 3)
}
