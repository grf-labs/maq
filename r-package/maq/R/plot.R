#' Plot the estimated Qini curve.
#'
#' Plot the estimated curve \eqn{Q(B), B \in (0, B_{max}]}. If the underlying estimated policy
#' \eqn{\pi_B} entails treating zero units, then this function returns an empty value.
#'
#' @param x A maq object.
#' @param ... Additional arguments passed to plot.
#' @param add Whether to add to an already existing plot. Default is FALSE.
#' @param horizontal.line Whether to draw a horizontal line where the Qini curve plateaus.
#'  Only applies if add = TRUE and the maq object is fit with a maximum `budget` that is sufficient
#'  to treat all units that are expected to benefit.
#'  Default is TRUE.
#' @param ci.args A list of optional arguments to `lines()` for drawing 95 % confidence bars.
#'  Set to NULL to ignore CIs.
#' @param grid.step The spend grid increment size to plot the curve on. Default is
#'  `max(floor(length(path.length) / 1000), 1)` where path.length is the size of the
#'  grid underlying the estimated Qini curve.
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
    if (x[["_path"]]$complete.path) {
      len <- length(spend)
      # Retrieve the main plot's end point (R by default extends the xrange by 4 percent)
      xmax <- graphics::par("usr")[2] / 1.04
      if (xmax > spend[len]) {
        spend <- c(spend, seq(spend[len], xmax, length.out = 100))
        gain <- c(gain, rep(gain[len], 100))
        std.err <- c(std.err, rep(std.err[len], 100))
      }
    }
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
