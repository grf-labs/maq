#' Title
#'
#' @param reward todo
#' @param cost todo
#' @param budget todo
#' @param R todo
#' @param sample.weights todo
#' @param clusters todo
#' @param tie.breaker todo
#' @param num.threads todo
#' @param seed todo
#'
#' @return todo
#'
#' @examples
#' \donttest{
#' # Train a
#' }
#'
#' @export
maq <- function(reward,
                cost,
                budget,
                R = 200,
                sample.weights = NULL,
                clusters = NULL,
                tie.breaker = NULL,
                num.threads = NULL,
                seed = runif(1, 0, .Machine$integer.max)) {
  if (NROW(reward) != NROW(cost) || NCOL(reward) != NCOL(cost)
        || anyNA(reward) || anyNA(cost)) {
    stop("reward and cost should be matrices of equal size with no missing values.")
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
    tie.breaker <- sample.int(NROW(reward))
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

  ret <- solver_rcpp(as.matrix(reward), as.matrix(cost), sample.weights, tie.breaker, clusters,
                     samples.per.cluster, budget, R, num.threads, seed)

  output <- list()
  class(output) <- "maq"
  output[["_path"]] <- ret
  output[["seed"]] <- seed
  output[["dim"]] <- c(NROW(cost), NCOL(cost))
  output[["budget"]] <- budget

  output
}

#' Get estimate of gain given a spend level.
#'
#' Gets estimates of
#'
#' @param object A maq object.
#' @param spend The spend level.
#'
#' @return An estimate of average gain along with standard errors.
#'
#' @examples
#' \donttest{
#' # Train a
#' }
#'
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

#' Predict optimal treatment allocation.
#'
#' Gets optimal alloction matrix for a given spend level.
#'
#' @param object A maq object.
#' @param spend The spend level.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A sparse matrix.
#'
#' @examples
#' \donttest{
#' # Train a
#' }
#'
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

#' Get
#'
#' Gets estimates of
#'
#' @param object A maq object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A thing
#'
#' @examples
#' \donttest{
#' # Train a
#' }
#'
#' @method summary maq
#' @export
summary.maq <- function(object,
                        ...) {

  data.frame(spend = object[["_path"]]$spend,
             gain = object[["_path"]]$gain,
             std.err = object[["_path"]]$std.err)
}

#' Print a maq object
#' @param x The maq object.
#' @param ... Additional arguments (currently ignored).
#'
#' @method print maq
#' @export
print.maq <- function(x, ...) {
  cat("MAQ object fit on", x$dim[1], "units and", x$dim[2], "arms with max budget", x$budget)
}
