#' Fit a Multi-Armed Qini.
#'
#' Consider \eqn{k = 1, \ldots, K} mutually exclusive and costly treatment arms,
#'  where k = 0 is a zero-cost control arm. Let \eqn{\hat \tau(\cdot)} be an _estimated_
#'  multi-armed treatment effect function and \eqn{C(\cdot)} a known cost function
#'  (where the k-th element of these vectors measures \eqn{E[Y_i(k) - Y_i(0) | X_i]} and
#'  \eqn{E[C_i(k) - C_i(0) | X_i]} where \eqn{Y_i(k)} are potential outcomes corresponding
#'  to the k-th treatment state, \eqn{C_i(k)} the cost of assigning unit i the k-th arm,
#'  and \eqn{X_i} a set of covariates). We provide estimates of the Qini curve:
#' \itemize{
#'    \item \eqn{Q(B) = E[\langle \pi_B(X_i), \tau(X_i)\rangle], B \in (0, B_{max}],}
#' }
#' which is the expected gain, at any budget constraint B, when assigning treatment in accordance
#'  to \eqn{\pi_B}, the treatment policy that optimally selects
#'  which arm to assign to which unit while incurring a cost less than or equal to B in expectation
#'  when using the given functions \eqn{\hat \tau(\cdot)} and \eqn{C(\cdot)}:
#' \itemize{
#'  \item \eqn{\pi_B = argmax_{\pi} \left\{E[\langle \pi(X_i), \hat \tau(X_i) \rangle]: E[\langle \pi(X_i), C(X_i) \rangle] \leq B \right\}.}
#' }
#' At a budget B, the k-th element of \eqn{\pi_B(X_i)} is 1 if assigning the k-th arm
#' to the i-th unit is optimal, and 0 otherwise.
#' The Qini curve can be used to quantify the value, as measured by the expected gain over
#'  assigning each unit the control arm when using the estimated function
#'  \eqn{\hat \tau(\cdot)} with cost structure \eqn{C(\cdot)} to allocate treatment,
#'  as we vary the available budget \eqn{B}.
#'
#'
#' @param reward A \eqn{n \cdot K} matrix of test set treatment effect estimates \eqn{\hat \tau(X_i)}.
#' (Note: the estimated function \eqn{\hat \tau(\cdot)} should be constructed on a held-out training set)
#' @param cost A \eqn{n \cdot K} matrix of test set costs \eqn{C(X_i) > 0}, where entry (i, k)
#'  measures the cost of assigning the i-th unit the k-th treatment arm.
#'  If the costs does not vary by unit, only by arm, this can also be a K-length vector.
#'  (Note: these costs need not be denominated on the same scale as the treatment effect estimates).
#' @param DR.scores An \eqn{n \cdot K} matrix of test set evaluation scores used to form an estimate of
#'  Q(B). With known treatment propensities \eqn{P[W_i|X_i]},
#'  these scores can be constructed via inverse-propensity weighting, i.e, with entry (i, k) equal to
#'  \eqn{\frac{\mathbf{1}(W_i=k)Y_i}{P[W_i=k | X_i]} - \frac{\mathbf{1}(W_i=0)Y_i}{P[W_i=0 | X_i]}}.
#'  In observational settings where \eqn{P[W_i|X_i]} has to be estimated, then an alternative is to
#'  construct these scores via augmented inverse-propensity weighting (AIPW) - yielding a doubly
#'  robust estimate of the Qini curve (for details, see the paper).
#' @param budget The maximum spend per unit, \eqn{B_{max}}, to fit the Qini curve on.
#'  Setting this to NULL (Default), will fit the path up to a maximum spend per unit
#'  where each unit that is expected to benefit (that is, \eqn{\hat \tau_k(X_i)>0}) is treated.
#' @param target.with.covariates If TRUE (Default), then the policy \eqn{\pi_B} takes covariates
#'  \eqn{X_i} into account. If FALSE, then the policy only takes the average reward
#'  \eqn{\bar \tau = E[\hat \tau(X_i)]} and average costs \eqn{\bar C = E[C(X_i)]} into account when
#'  allocating treatment. This can be used to construct a baseline Qini curve to assess the value
#'  of treatment targeting based on covariates.
#' @param R Number of bootstrap replicates for computing standard errors. Default is 0
#'  (only point estimates are computed).
#' @param paired.inference Whether to allow for paired tests with other Qini curves fit on the same
#'  evaluation data. If TRUE (Default) then the path of bootstrap replicates are stored in order to perform
#'  paired comparisons that account for the correlation between curves evaluated on the same data. This
#'  takes memory on the order of \eqn{O(RnK)} and requires the comparison objects to be fit with
#'  the same seed and R values as well as the same number of samples.
#' @param sample.weights Weights given to an observation in estimation.
#'  If NULL, each observation is given the same weight. Default is NULL.
#' @param clusters Vector of integers or factors specifying which cluster each observation corresponds to,
#'  which are used to construct clustered standard errors.
#'  Default is NULL (ignored).
#' @param tie.breaker An optional permutation of the integers 1 to n used to
#'  break potential ties in the optimal treatment allocation
#'  (only relevant if \eqn{\hat \tau(X)} takes on the same values for different samples
#'  \eqn{X_i} and \eqn{X_j}).
#'  If NULL, the ties are broken by the lowest sample id (i.e. the sample appearing first in the data).
#'  Default is NULL.
#' @param num.threads Number of threads used in bootstrap replicates. By default, the number of threads
#'  is set to the maximum hardware concurrency.
#' @param seed The seed of the C++ random number generator. Default is 42.
#'
#' @return A fit maq object.
#'
#' @references Sverdrup, Erik, Han Wu, Susan Athey, and Stefan Wager.
#'  "Qini Curves for Multi-Armed Treatment Rules".
#'  arXiv preprint arXiv:2306.11979, 2023.
#'
#' @examples
#' \donttest{
#' if (require("grf", quietly = TRUE)) {
#' # Fit a CATE estimator on a training sample.
#' n <- 3000
#' p <- 5
#' X <- matrix(runif(n * p), n, p)
#' W <- as.factor(sample(c("0", "1", "2"), n, replace = TRUE))
#' Y <- X[, 1] + X[, 2] * (W == "1") + 1.5 * X[, 3] * (W == "2") + rnorm(n)
#' train <- sample(1:n, n/2)
#'
#' tau.forest <- grf::multi_arm_causal_forest(X[train, ], Y[train], W[train])
#'
#' # Predict CATEs on held out evaluation data.
#' test <- -train
#' tau.hat <- predict(tau.forest, X[test, ], drop = TRUE)$predictions
#'
#' # Assume costs equal a unit's pre-treatment covariate - the following are a toy example.
#' cost <- cbind(X[test, 4] / 4, X[test, 5])
#'
#' # Fit an evaluation forest to compute doubly robust scores on the test set.
#' eval.forest <- grf::multi_arm_causal_forest(X[test, ], Y[test], W[test])
#' DR.scores <- grf::get_scores(eval.forest, drop = TRUE)
#'
#' # Fit a Qini curve on evaluation data, using 200 bootstrap replicates for confidence intervals.
#' ma.qini <- maq(tau.hat, cost, DR.scores, R = 200)
#'
#' # Plot the Qini curve.
#' plot(ma.qini)
#' legend("topleft", c("All arms", "95% CI"), lty = c(1, 3))
#'
#' # Get an estimate of gain at a given spend per unit along with standard errors.
#' average_gain(ma.qini, spend = 0.2)
#'
#' # Get the treatment allocation matrix at a given spend per unit.
#' pi.mat <- predict(ma.qini, spend = 0.2)
#'
#' # If the treatment randomization probabilities are known, then an alternative to
#' # evaluation via AIPW scores is to use inverse-propensity weighting (IPW).
#' W.hat <- rep(1/3, 3)
#' IPW.scores <- get_ipw_scores(Y[test], W[test], W.hat)
#' mq.ipw <- maq(tau.hat, cost, IPW.scores)
#'
#' plot(mq.ipw, add = TRUE, col = 2)
#' legend("topleft", c("All arms", "95% CI", "All arms (IPW)"), col = c(1, 1, 2), lty = c(1, 3, 1))
#'
#' # Estimate some baseline policies.
#' # a) A policy that ignores covariates and only takes the average reward/cost into account.
#' qini.avg <- maq(tau.hat, cost, DR.scores, target.with.covariates = FALSE, R = 200)
#'
#' # b) A policy that only use arm 1.
#' qini.arm1 <- maq(tau.hat[, 1], cost[, 1], DR.scores[, 1], R = 200)
#'
#' # c) A policy that only use arm 2.
#' qini.arm2 <- maq(tau.hat[, 2], cost[, 2], DR.scores[, 2], R = 200)
#'
#' plot(ma.qini, ci.args = NULL)
#' plot(qini.avg, col = 2, add = TRUE, ci.args = NULL)
#' plot(qini.arm1, col = 3, add = TRUE, ci.args = NULL)
#' plot(qini.arm2, col = 4, add = TRUE, ci.args = NULL)
#' legend("topleft", c("All arms (targeting)", "All arms (without targeting)", "Arm 1", "Arm 2"),
#'        col = 1:4, lty = 1)
#'
#' # Estimate the value of employing all arms over a random allocation.
#' difference_gain(ma.qini, qini.avg, spend = 0.2)
#'
#' # Estimate the value of targeting with both arms as opposed to targeting with only arm 1.
#' difference_gain(ma.qini, qini.arm1, spend = 0.2)
#'
#' # Estimate the value of targeting with both arms as opposed to targeting with only arm 2.
#' difference_gain(ma.qini, qini.arm2, spend = 0.2)
#'
#' # Compare treatment strategies over a range of budget values by estimating an area between
#' # two curves up to a given spend point.
#' integrated_difference(ma.qini, qini.arm1, spend = 0.3)
#' }
#' }
#'
#' @export
maq <- function(reward,
                cost,
                DR.scores,
                budget = NULL,
                target.with.covariates = TRUE,
                R = 0,
                paired.inference = TRUE,
                sample.weights = NULL,
                clusters = NULL,
                tie.breaker = NULL,
                num.threads = NULL,
                seed = 42) {
  if (NROW(reward) != NROW(DR.scores) || NCOL(reward) != NCOL(DR.scores)
        || anyNA(reward) || anyNA(cost) || anyNA(DR.scores)) {
    stop("reward, costs, and evaluation scores should have conformable dimension, with no missing values.")
  }
  if (is.vector(cost) && length(cost) == NCOL(reward)) {
    cost <- matrix(cost, 1, length(cost), byrow = TRUE)
  } else {
    if (NROW(cost) != NROW(reward) || NCOL(cost) != NCOL(reward)) {
      stop("reward, costs, and evaluation scores should have conformable dimension, with no missing values.")
    }
  }

  if (any(cost <= 0)) {
    stop("Costs should be > 0.")
  }

  if (is.null(budget)) {
    max.budget <- .Machine$double.xmax
  } else {
    max.budget <- budget
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

  ret <- solver_rcpp(as.matrix(reward), as.matrix(DR.scores), as.matrix(cost),
                     sample.weights, tie.breaker, clusters,
                     max.budget, target.with.covariates, paired.inference, R, num.threads, seed)

  output <- list()
  class(output) <- "maq"
  output[["_path"]] <- ret
  output[["seed"]] <- seed
  output[["target.with.covariates"]] <- target.with.covariates
  output[["paired.inference"]] <- paired.inference
  output[["R"]] <- R
  output[["dim"]] <- c(NROW(reward), NCOL(reward))
  output[["budget"]] <- if (is.null(budget)) {
    max(ret$spend[length(ret$spend)], 0)
  } else {
    budget
  }

  output
}

#' Predict treatment allocation.
#'
#' Get an estimate of the policy \eqn{\pi_B(X_i)} at a spend level B.
#'  \eqn{\pi_B(X_i)} is a K-dimensional vector where the k-th element is 1 if assigning the k-th
#'  arm to unit i is optimal at a given spend B, and 0 otherwise (with all entries 0 if the
#'  control arm is assigned).
#'  Depending on the value of B, \eqn{\pi_B(X_j)} might be fractional for at most one unit j.
#'  There are two such cases - the first one is when there is not sufficient budget left to assign j an
#'  initial arm. The second is if there is not sufficient budget to upgrade unit j from arm k to k'.
#'  In these cases \eqn{\pi_B(X_j)} takes on one, or two fractional values, respectively,
#'  representing an assignment probability of a given arm.
#'
#'
#' @param object A maq object.
#' @param spend The spend level B.
#' @param type If "matrix" (Default), then return a matrix where the i-th entry equals
#'  \eqn{\pi_B(X_i)} as described above.
#'  If "vector", then \eqn{\pi_B(X_i)} is instead encoded taking values in the set \{0, 1, ..., K\}.
#'  If the allocation is fractional at the given B, this option returns the policy corresponding
#'  to the previous/lower value of the spend path, at which point the policy is integer-valued, but
#'  incurs a cost less than B in expectation.
#'
#' @param ... Additional arguments (currently ignored).
#'
#' @return A matrix with row i equal to \eqn{\pi_B(X_i)}. If `type = "vector"` then an
#'  n-length vector with elements equal to the arm (from 0 to K) that is assigned at the given spend B
#'  (note: if the treatment allocation contains a fractional entry at the given B, then the returned
#'  vector is the policy at the nearest spend B' in the solution path where the allocation is
#'  integer-valued but incurs a cost B' < B).
#'
#' @examples
#' \donttest{
#' # Generate some toy data and fit a solution path.
#' n <- 10
#' K <- 4
#' reward <- matrix(rnorm(n * K), n, K)
#' cost <- matrix(runif(n * K), n, K)
#' DR.scores <- reward + rnorm(n)
#' path <- maq(reward, cost, DR.scores)
#'
#' # Get the treatment allocation matrix
#' pi.mat <- predict(path, 0.1)
#' pi.mat
#' # pi.mat might have fractional entries for a single unit but satisfies
#' # the budget in expectation exactly.
#' sum(cost * pi.mat) / n
#'
#' # Get the treatment allocation instead encoded in the set \{0, 1, ..., K\}.
#' pi.vec <- predict(path, 0.1, type = "vector")
#' pi.vec
#' # If a unit has a fractional entry, then pi.vec will incur a cost slightly
#' # lower than 0.1.
#' sum(cost[cbind(1:n, pi.vec)]) / n
#'
#' # Retrieve the underlying solution path.
#' data.path <- summary(path)
#' # If we predict at a spend level on this grid, say entry 5,
#' # then the policy is integer-valued:
#' spend <- data.path$spend[5]
#' predict(path, spend)
#' predict(path, spend, type = "vector")
#' }
#'
#' @method predict maq
#' @export
predict.maq <- function(object,
                        spend,
                        type = c("matrix", "vector"),
                        ...) {
  type <- match.arg(type)
  if (!object[["_path"]]$complete.path && spend > object$budget) {
    stop("maq path is not fit beyond given spend level.")
  }

  spend.grid <- object[["_path"]]$spend
  path.idx <- findInterval(spend, spend.grid) # nearest path index (lower bound)
  if (path.idx == 0) {
    if (type == "matrix") {
      return (matrix(0, object[["dim"]][1], object[["dim"]][2]))
    } else {
      return (rep(0, object[["dim"]][1]))
    }
  }

  ipath <- object[["_path"]]$ipath[1:path.idx] + 1 # +1: R index.
  kpath <- object[["_path"]]$kpath[1:path.idx] + 1
  ix <- !duplicated(ipath, fromLast = TRUE)

  if (type == "matrix") {
    pi.mat <- matrix(0, object[["dim"]][1], object[["dim"]][2])
    pi.mat[cbind(ipath[ix], kpath[ix])] <- 1
  } else {
    pi.vec <- rep(0, object[["dim"]][1])
    pi.vec[ipath[ix]] <- kpath[ix]
    return (pi.vec)
  }

  if (path.idx == length(spend.grid)) {
    return (pi.mat)
  }

  # fractional adjustment? (only done when return type is a matrix)
  spend.diff <- spend - spend.grid[path.idx]
  fraction <- spend.diff / (spend.grid[path.idx + 1] - spend.grid[path.idx])

  next.unit <- object[["_path"]]$ipath[path.idx + 1] + 1
  next.arm <- object[["_path"]]$kpath[path.idx + 1] + 1
  prev.arm <- which(pi.mat[next.unit, ] == 1) # already assigned an arm?

  pi.mat[next.unit, next.arm] <- fraction
  if (length(prev.arm) > 0) {
    pi.mat[next.unit, prev.arm] <- 1 - fraction
  }

  pi.mat
}

#' Get estimate of gain given a spend level.
#'
#' Get an estimate of Q(B).
#'
#' @param object A maq object.
#' @param spend The spend level B.
#'
#' @return An estimate of Q(B) along with standard errors.
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

#' Get estimate of difference in gain given a spend level with paired standard errors.
#'
#' Given two Qini curves, \eqn{Q_a} and \eqn{Q_b}, get an estimate of the difference
#' \eqn{Q_a(B) - Q_b(B)}, at a spend level B.
#'
#' @param object.lhs A maq object \eqn{Q_a} to subtract from.
#' @param object.rhs A maq object \eqn{Q_b} to subtract with.
#' @param spend The spend level B.
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

  point.estimate <- average_gain(object.lhs, spend)[[1]] - average_gain(object.rhs, spend)[[1]]
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

  c(estimate = point.estimate, std.err = std.err)
}

#' Get estimate of area the between two Qini curves with paired standard errors.
#'
#' Given two Qini curves, \eqn{Q_a} and \eqn{Q_b}, and a maximum spend \eqn{\bar B},
#'  get an estimate of the integrated difference
#'  \eqn{\int_{0}^{\bar B} (Q_a(B) - Q_b(B))dB}.
#'
#' @param object.lhs A maq object \eqn{Q_a} to subtract from.
#' @param object.rhs A maq object \eqn{Q_b} to subtract with.
#' @param spend The spend level \eqn{\bar B}.
#'
#' @return An estimate of the area between the two curves along with standard errors.
#' @export
integrated_difference <- function(object.lhs,
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
  R <- object.lhs[["R"]]

  # Estimate an AUC via estimating the difference \int_{0}^{\bar B} Q_a(B)dB - \int_{0}^{\bar B} Q_b(B)dB.
  .get_estimate <- function(gain.path, spend.grid, path.idx) {
    if (path.idx == 0) {
      estimate <- 0
    } else if (abs(spend.grid[length(spend.grid)] - spend) < 1e-10 || path.idx == length(spend.grid)) {
      area.offset <- 0
      # Are we summing beyond the point at which the curve plateaus?
      if (spend > spend.grid[path.idx]) {
        spend.delta <- spend - spend.grid[path.idx]
        area.offset <- spend.grid[path.idx] * spend.delta
      }
      estimate <- mean(gain.path, na.rm = TRUE) + area.offset
    } else {
      interp.ratio <- (spend - spend.grid[path.idx]) / (spend.grid[path.idx + 1] - spend.grid[path.idx])
      estimate <- mean(gain.path[1:path.idx], na.rm = TRUE) +
       (gain.path[path.idx + 1] - gain.path[path.idx]) * interp.ratio
    }

    estimate
  }

  gain.path.lhs <- object.lhs[["_path"]]$gain
  spend.grid.lhs <- object.lhs[["_path"]]$spend
  path.idx.lhs <- findInterval(spend, spend.grid.lhs)

  gain.path.rhs <- object.rhs[["_path"]]$gain
  spend.grid.rhs <- object.rhs[["_path"]]$spend
  path.idx.rhs <- findInterval(spend, spend.grid.rhs)

  point.estimate <- .get_estimate(gain.path.lhs, spend.grid.lhs, path.idx.lhs) -
    .get_estimate(gain.path.rhs, spend.grid.rhs, path.idx.rhs)

  # Compute paired std.errors
  estimates.lhs <- unlist(lapply(seq_len(R), function(bs) {
    .get_estimate(object.lhs[["_path"]]$gain.bs[[bs]], spend.grid.lhs, path.idx.lhs)
  }))
  estimates.rhs <- unlist(lapply(seq_len(R), function(bs) {
    .get_estimate(object.rhs[["_path"]]$gain.bs[[bs]], spend.grid.rhs, path.idx.rhs)
  }))
  std.err <- stats::sd(estimates.lhs - estimates.rhs, na.rm = TRUE)
  if (is.na(std.err)) {
    std.err <- 0
  }

  c(estimate = point.estimate, std.err = std.err)
}
