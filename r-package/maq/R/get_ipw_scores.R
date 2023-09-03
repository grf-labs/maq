#' Construct evaluation scores via inverse-propensity weighting.
#'
#' Construct an evaluation score matrix via IPW, where entry (i, k) equals
#' \itemize{
#'  \item \eqn{\frac{\mathbf{1}(W_i=k)Y_i}{P[W_i=k | X_i]} - \frac{\mathbf{1}(W_i=0)Y_i}{P[W_i=0 | X_i]}},
#' }
#' where \eqn{W_i} is the treatment assignment of unit i and \eqn{Y_i} the observed outcome.
#' \eqn{k = 1 \ldots K} are one of K treatment arms and k = 0 is the control arm.
#'
#'
#' @param Y The observed outcome.
#' @param W The observed treatment assignment (must be a factor vector,
#'  where the first factor level is the control arm).
#' @param W.hat Optional treatment propensities. If these vary by unit and arm, then
#'  this should be a matrix where with the treatment assignment
#'  probability of units to arms, with columns corresponding to the levels of `W`.
#'  If these only vary by arm, a vector can also be supplied.
#'  If W.hat is NULL (Default), then the assignment probabilities are assumed to be uniform
#'  and the same for each arm.
#'
#' @return An \eqn{n \cdot K} matrix of evaluation scores.
#'
#' @examples
#' \donttest{
#'
#' # Draw some equally likely samples from control arm A and treatment arms B and C.
#' n <- 5000
#' W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE))
#' Y <- 42 * (W == "B") - 42 * (W == "C") + rnorm(n)
#' IPW.scores <- get_ipw_scores(Y, W)
#' # Should be ~ 42 and -42.
#' colMeans(IPW.scores)
#'
#' # Draw non-uniformly from the different arms.
#' W.hat <- c(0.2, 0.2, 0.6)
#' W <- as.factor(sample(c("A", "B", "C"), n, replace = TRUE, prob = W.hat))
#' Y <- 42 * (W == "B") - 42 * (W == "C") + rnorm(n)
#' IPW.scores <- get_ipw_scores(Y, W, W.hat = W.hat)
#' # Should still be ~ 42 and -42.
#' colMeans(IPW.scores)
#' }
#'
#' @export
get_ipw_scores <- function(Y,
                           W,
                           W.hat = NULL) {
  if (!is.factor(W)) {
    stop("W should be a factor vector.")
  }
  if (length(Y) != length(W) || anyNA(Y) || anyNA(W)) {
    stop("Y and W should be equal-length vectors with no missing entries.")
  }
  K.plus1 <- nlevels(W)
  if (is.null(W.hat)) {
    W.hat <- rep(1 / K.plus1, K.plus1)
  }
  if (anyNA(W.hat) || any(W.hat <= 0) || any(W.hat >= 1)) {
    stop("W.hat entries should be non-missing and between (0, 1).")
  }
  if (is.vector(W.hat) && length(W.hat) == K.plus1) {
    W.hat <- matrix(W.hat, length(Y), length(W.hat), byrow = TRUE)
  }
  if (NROW(W.hat) != length(W) || NCOL(W.hat) != K.plus1) {
    stop("W.hat should either be a (K+1)-length vector or a n*(K+1) matrix of treatment propensities.")
  }
  if (any(abs(rowSums(W.hat) - 1) > 1e-5)) {
    stop("W.hat propensities should sum to 1.")
  }

  observed.W <- match(W, levels(W))
  Y.mat <- matrix(0, length(W), K.plus1)
  Y.mat[cbind(seq_along(observed.W), observed.W)] <- Y
  Y.ipw <- Y.mat / W.hat

  Y.ipw[, -1] - Y.ipw[, 1]
}
