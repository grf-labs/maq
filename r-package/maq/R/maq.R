#' Title
#'
#' @param reward todo
#' @param cost todo
#' @param budget todo
#' @param R todo
#' @param sample.weights todo
#' @param clusters todo
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
                R = 150,
                sample.weights = NULL,
                clusters = NULL,
                num.threads = NULL,
                seed = runif(1, 0, .Machine$integer.max)) {
  if (NROW(reward) != NROW(cost) || NCOL(reward) != NCOL(cost)
        || anyNA(reward) || anyNA(cost) || NCOL(reward) < 2) {
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
  } else if (length(sample.weights) != NROW(reward) || anyNA(sample.weights)) {
    stop("sample.weights have incorrect length.")
  }
  # todo add wieights > 0 and say drop instead of giving wt =0...

  if (is.null(clusters) || length(clusters) == 0) {
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
  if (is.null(num.threads)) {
    num.threads <- 0
  } else if (num.threads < 0) {
    stop("num.threads should be a non-negative integer.")
  }
  if (!is.numeric(seed) || seed < 0) {
    stop("seed should be a non-negative integer.")
  }

  # create a matrix of sort indices ordered by increasing cost
  # the following is just a faster way of sorting K items n times, like
  # `ix.order.slow <- t(apply(cost, 1, order))`
  # by instead doing a grouped sort of size n * K.
  # ix.order <- matrix(order(col(t(cost)), t(cost), -1 * t(reward)), ncol = ncol(cost), byrow = TRUE) -
    # ncol(cost) * (row(cost) - 1)
  ix.order <- matrix(order(col(t(cost)), t(cost)), ncol = ncol(cost), byrow = TRUE) -
    ncol(cost) * (row(cost) - 1) - 1 # -1: C++ index

  ret <- solver_rcpp(as.matrix(reward), as.matrix(cost), ix.order, sample.weights, clusters,
                     samples.per.cluster, budget, R, num.threads, seed)

# browser()
# print(system.time({
  if (length(ret[["t0"]]$spend) > 5) { # TODO
    t.grid <- lapply(seq_len(R), function(i) {
      approx(ret[["t"]]$spend[[i]],
             ret[["t"]]$gain[[i]],
             ret[["t0"]]$spend,
             rule = 2, #should really be gain[min]... on left?
             ties = "ordered")$y
    })
  } else {
    t.grid <- list()
  }
  if (length(t.grid) > 0) {
    t.mat <- matrix(unlist(t.grid), R, length(t.grid[[1]]), byrow = TRUE)
    if (requireNamespace("matrixStats", quietly = TRUE)) {
      std.err <- matrixStats::colSds(t.mat) # TODO eps
    } else {
      std.err <- apply(t.mat, 2, sd)
    }
  } else {
    std.err <- 0
  }
# }))

  output <- list()
  class(output) <- "maq"
  output[["_path"]] <- ret[["t0"]]
  output[["_path"]][["std.err"]] <- std.err
  output[["_path.bs"]] <- ret[["t"]]
  output[["seed"]] <- seed
  output[["dim"]] <- c(NROW(cost), NCOL(cost))
  output[["budget"]] <- budget

  output
}

#' Predict with a
#'
#' Gets optimal alloc matrix for a given spend level.
#'
#' @param object The trained
#' @param ... Additional arguments (currently ignored).
#'
#' @return Vector of
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
  ## nearest path index
  path.idx <- findInterval(spend, object[["_path"]]$spend)
  if (path.idx == 0) {
    path.idx <- 1
  }
  ipath <- object[["_path"]]$ipath[1:path.idx] + 1 # +1: R index.
  kpath <- object[["_path"]]$kpath[1:path.idx] + 1
  ix <- !duplicated(ipath, fromLast = TRUE)

  Matrix::sparseMatrix(ipath[ix], kpath[ix], dims = object[["dim"]])
}

#' Get estimates of
#'
#' Gets es
#'
#' @param object The traine
#' @param ... Additional arguments (currently ignored).
#'
#' @return A thing
#'
#' @examples
#' \donttest{
#' # Train a
#' }
#'
#' @method coef maq
#' @export
coef.maq <- function(object,
                     blabl,
                     ...) {
  ##
  cbind(spend = object[["_path"]]$spend, gain = object[["_path"]]$gain)
}

#' Print a maq object
#' @param x The
#' @param ... Additional arguments (currently ignored).
#'
#' @method print maq
#' @export
print.maq <- function(x, ...) {
  cat("hey")
}

# summary
