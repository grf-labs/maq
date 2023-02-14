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
  if (is.null(num.threads)) {
    num.threads <- 0
  } else if (num.threads < 0) {
    stop("num.threads should be a non-negative integer.")
  }
  if (!is.numeric(seed) || seed < 0) {
    stop("seed should be a non-negative integer.")
  }

  # TODO /n at the end istead avoid 2 *O(n*k)?
  ret <- solver_rcpp(as.matrix(reward) / NROW(reward), as.matrix(cost) / NROW(cost), sample.weights, clusters,
                     samples.per.cluster, budget, R, num.threads, seed)

  output <- list()
  class(output) <- "maq"
  output[["_path"]] <- ret
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
                        type = c("gain", "pi.matrix"),
                        ...) {
  type <- match.arg(type)
  ## nearest path index
  spend.grid <- object[["_path"]]$spend
  path.idx <- findInterval(spend, spend.grid)

  if (type == "gain") {
    gain.path <- object[["_path"]]$gain
    se.path <- object[["_path"]]$std.err
    if (path.idx == 0) {
      estimate <- NA
      std.err <- NA
    } else if (path.idx == length(spend.grid)) {
      estimate <- spend.grid[path.idx]
      std.err <- se.path[path.idx]
    } else {
      interp.ratio <- (spend - spend.grid[path.idx]) / (spend.grid[path.idx + 1] - spend.grid[path.idx])
      estimate <- gain.path[path.idx] + (gain.path[path.idx + 1] - gain.path[path.idx]) * interp.ratio
      std.err <- se.path[path.idx] + (se.path[path.idx + 1] - se.path[path.idx]) * interp.ratio
    }

    return (c(estimate = estimate, std.err = std.err))
  }

  if (path.idx == 0) {
    path.idx <- 1 # TODO
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
