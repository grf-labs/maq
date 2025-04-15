#' Plot the estimated Qini curve.
#'
#' Plot the estimated curve \eqn{Q(B), B \in (0, B_{max}]}. If the underlying estimated policy
#' \eqn{\pi_B} entails treating zero units (that is, all the estimated treatment effects are
#'  negative) then this function returns an empty value.
#'
#' @param x A maq object.
#' @param ... Additional arguments passed to plot.
#' @param add Whether to add to an already existing plot. Default is FALSE.
#' @param horizontal.line Whether to draw a horizontal line where the Qini curve plateaus.
#'  Only applies if the maq object is fit with a maximum `budget` that is sufficient
#'  to treat all units that are expected to benefit.
#'  Default is TRUE.
#' @param ci.args A list of optional arguments to `lines()` for drawing 95 % confidence bars.
#'  Set to NULL to ignore CIs.
#' @param grid.step The spend grid increment size to plot the curve on. Default is
#'  `max(floor(length(path.length) / 1000), 1)` where path.length is the size of the
#'  grid underlying the estimated Qini curve.
#'
#' @return A data.frame with the data making up the plot (point estimates and lower/upper 95% CIs)
#'
#' @examples
#' \donttest{
#' # Generate toy data and customize plots.
#' n <- 500
#' K <- 1
#' reward <- matrix(1 + rnorm(n * K), n, K)
#' scores <- reward + matrix(rnorm(n * K), n, K)
#' cost <- 1
#'
#' # Fit Qini curves.
#' qini <- maq(reward, cost, scores, R = 200)
#' qini.avg <- maq(reward, cost, scores, R = 200, target.with.covariates = FALSE)
#'
#' # The plot method invisibly returns the plot data as a data.frame,
#' # which allows for custom plotting with external libraries.
#' df.qini.baseline <- plot(qini.avg)
#' df.qini <- plot(qini, add = TRUE, col = "red")
#'
#' # Make an alternate plot style, using, for example, ggplot.
#' if (require("ggplot2", quietly = TRUE)) {
#' ggplot(df.qini, aes(x = spend, y = gain)) +
#'   geom_ribbon(aes(ymin = gain - 1.96 * std.err,
#'                   ymax = gain + 1.96 * std.err),
#'               fill = "lightgray") +
#'   geom_line(linewidth = 2) +
#'   ylab("Policy value") +
#'   xlab("Fraction treated") +
#'   geom_line(data = df.qini.baseline, aes(x = spend, y = gain), lty = 2)
#' }
#'
#' # `scale_maq()` rescales policy gain curves for specific application units.
#' # Plot policy values for a maximum allocation of, for example, 500 units.
#' plot(scale_maq(qini, 500), xlab = "Units treated")
#' }
#' @method plot maq
#' @export
plot.maq <- function(x,
                     ...,
                     add = FALSE,
                     horizontal.line = TRUE,
                     ci.args = list(),
                     grid.step = NULL
                     ) {
  new.args <- list(...)
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
  if (horizontal.line && x[["_path"]]$complete.path) {
    len <- length(spend)
    # Are we creating a new plot with a user-supplied xlim that extends beyond the fit curve?
    if (!add && "xlim" %in% names(new.args)) {
      xmax <- new.args$xlim[2]
      if (xmax > spend[len]) {
        spend <- c(spend, seq(spend[len], xmax, length.out = 100))
        gain <- c(gain, rep(gain[len], 100))
        std.err <- c(std.err, rep(std.err[len], 100))
        }
    }
    # Else, we're adding to an existing plot
    if (add) {
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
  plot.args[names(new.args)] <- new.args

  lines.args <- list(lty = 3, col = plot.args$col)
  lines.args[names(ci.args)] <- ci.args

  if (!add || grDevices::dev.cur() == 1L) {
    do.call(plot, c(list(x = spend, y = gain), plot.args))
  } else {
    do.call(graphics::lines, c(list(x = spend, y = gain), plot.args))
  }

  if (!is.null(ci.args) && x[["R"]] > 1) {
    do.call(graphics::lines, c(list(x = spend, y = lb), lines.args))
    do.call(graphics::lines, c(list(x = spend, y = ub), lines.args))
  }

  invisible(data.frame(spend, gain, std.err, lb, ub))
}
