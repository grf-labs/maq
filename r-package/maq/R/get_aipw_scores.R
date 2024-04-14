#' Construct evaluation scores via augmented inverse-propensity weighting.
#'
#' A simple convenience function to construct an AIPW-based evaluation score given estimates
#'  of conditional means and treatment propensities.
#'
#' @param Y The observed outcome.
#' @param W The observed treatment assignment (must be a factor vector,
#'  where the first factor level is the control arm).
#' @param Y.hat A matrix of conditional mean estimates for each arm, \eqn{E[Y_i | W_i = k, X_i]}.
#' @param W.hat Optional treatment propensities. If these vary by unit and arm, then
#'  this should be a matrix with the treatment assignment
#'  probability of units to arms, with columns corresponding to the levels of `W`.
#'  If these only vary by arm, a vector can also be supplied.
#'  If W.hat is NULL (Default), then the assignment probabilities are assumed to be uniform
#'  and the same for each arm.
#'
#' @return An \eqn{n \cdot K} matrix of evaluation scores
#'  (eqn (13) in the multi-armed Qini paper).
#'
#' @references Robins, James M, Andrea Rotnitzky, and Lue Ping Zhao.
#'  "Estimation of regression coefficients when some regressors are not always observed."
#'  Journal of the American statistical Association, 89(427), 1994.
#'
#' @references Sverdrup, Erik, Han Wu, Susan Athey, and Stefan Wager.
#'  "Qini Curves for Multi-Armed Treatment Rules".
#'  arXiv preprint arXiv:2306.11979, 2023.
#'
#' @examples
#' \donttest{
#' if (require("grf", quietly = TRUE)) {
#' n <- 3000
#' p <- 5
#' X <- matrix(runif(n * p), n, p)
#' W <- as.factor(sample(c("0", "1", "2"), n, replace = TRUE))
#' Y <- X[, 1] + X[, 2] * (W == "1") + 1.5 * X[, 3] * (W == "2") + rnorm(n)
#'
#' # Fit a CATE estimator on a training sample.
#' train <- sample(1:n, n/2)
#' tau.forest <- grf::multi_arm_causal_forest(X[train, ], Y[train], W[train])
#'
#' # Predict CATEs on held out evaluation data.
#' test <- -train
#' tau.hat <- predict(tau.forest, X[test, ], drop = TRUE)$predictions
#' # Form costs.
#' cost <- cbind(X[test, 4] / 4, X[test, 5])
#'
#' # Estimate nuisance components for test set AIPW scores.
#' X.test <- X[test, ]
#' Y.test <- Y[test]
#' W.test <- W[test]
#'
#' # Fit models for E[Y | W = k, X], k = 0, 1, 2, using for example separate random forests.
#' Y0.forest <- grf::regression_forest(X.test[W.test == 0, ], Y.test[W.test == 0])
#' Y1.forest <- grf::regression_forest(X.test[W.test == 1, ], Y.test[W.test == 1])
#' Y2.forest <- grf::regression_forest(X.test[W.test == 2, ], Y.test[W.test == 2])
#' Y.hat = cbind(
#'    mu0 = predict(Y0.forest, X.test)$predictions,
#'    mu1 = predict(Y1.forest, X.test)$predictions,
#'    mu2 = predict(Y2.forest, X.test)$predictions
#' )
#'
#' # If unknown, estimate the propensity scores E[W = k | X].
#' W.hat <- predict(grf::probability_forest(X.test, W.test))$predictions
#'
#' # Form doubly robust scores.
#' DR.scores <- get_aipw_scores(Y.test, W.test, Y.hat, W.hat)
#'
#' # Fit a Qini curve estimated with forest-based AIPW.
#' qini <- maq(tau.hat, cost, DR.scores, R = 200)
#' plot(qini)
#' }
#' }
#'
#' @export
get_aipw_scores <- function(Y,
                            W,
                            Y.hat,
                            W.hat = NULL) {
  if (!is.factor(W)) {
    stop("W should be a factor vector.")
  }
  if (length(Y) != length(W) || anyNA(Y) || anyNA(W)) {
    stop("Y and W should be equal-length vectors with no missing entries.")
  }
  K.plus1 <- nlevels(W)
  if (K.plus1 == 1) {
    stop("There is only one treatment level.")
  }
  if (is.null(W.hat)) {
    W.hat <- rep(1 / K.plus1, K.plus1)
  }
  if (is.vector(W.hat) && length(W.hat) == K.plus1) {
    W.hat <- matrix(W.hat, length(W), length(W.hat), byrow = TRUE)
  }
  if (NROW(W.hat) != length(W) || NCOL(W.hat) != K.plus1) {
    stop("W.hat should either be a (K+1)-length vector or a n*(K+1) matrix of treatment propensities.")
  }
  if (anyNA(W.hat) || any(W.hat < 0) || any(W.hat > 1)) {
    stop("W.hat entries should be non-missing and between (0, 1).")
  }
  if (any(abs(rowSums(W.hat) - 1) > 1e-5)) {
    stop("W.hat propensities should sum to 1.")
  }
  if (!identical(dim(Y.hat), dim(W.hat)) || anyNA(Y.hat)) {
    stop("Y.hat should be a matrix of conditional mean estimates for each arm.")
  }
  if (any(W.hat < 0.05) || any(W.hat > 0.95)) {
    warning("Some treatment propensities are lower/higher than 0.05/0.95 - overlap may be an issue.")
  }

  observed.W <- match(W, levels(W))
  observed.W.idx <- cbind(seq_along(W), observed.W)

  IPW <- matrix(0, length(W), K.plus1)
  IPW[observed.W.idx] <- 1 / W.hat[observed.W.idx]
  control <- IPW[, 1] > 0
  IPW[control, -1] <- -1 * IPW[control, 1]
  IPW <- IPW[, -1, drop = FALSE]

  out <- Y.hat[, -1, drop = FALSE] - Y.hat[, 1] + (Y - Y.hat[observed.W.idx]) * IPW
  colnames(out) <- paste(levels(W)[-1], "-", levels(W)[1])

  out
}
