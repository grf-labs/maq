% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.maq}
\alias{plot.maq}
\title{Plot the estimated Qini curve.}
\usage{
\method{plot}{maq}(
  x,
  ...,
  add = FALSE,
  horizontal.line = TRUE,
  ci.args = list(),
  grid.step = NULL
)
}
\arguments{
\item{x}{A maq object.}

\item{...}{Additional arguments passed to plot.}

\item{add}{Whether to add to an already existing plot. Default is FALSE.}

\item{horizontal.line}{Whether to draw a horizontal line where the Qini curve plateaus.
Only applies if the maq object is fit with a maximum \code{budget} that is sufficient
to treat all units that are expected to benefit.
Default is TRUE.}

\item{ci.args}{A list of optional arguments to \code{lines()} for drawing 95 \% confidence bars.
Set to NULL to ignore CIs.}

\item{grid.step}{The spend grid increment size to plot the curve on. Default is
\code{max(floor(length(path.length) / 1000), 1)} where path.length is the size of the
grid underlying the estimated Qini curve.}
}
\value{
A data.frame with the data making up the plot (point estimates and lower/upper 95\% CIs)
}
\description{
Plot the estimated curve \eqn{Q(B), B \in (0, B_{max}]}. If the underlying estimated policy
\eqn{\pi_B} entails treating zero units (that is, all the estimated treatment effects are
negative) then this function returns an empty value.
}
\examples{
\donttest{
if (require("ggplot2", quietly = TRUE)) {
# Generate toy data and customize plots.
n <- 500
K <- 1
reward <- matrix(1 + rnorm(n * K), n, K)
scores <- reward + matrix(rnorm(n * K), n, K)
cost <- 1

# Fit Qini curves.
qini.avg <- maq(reward, cost, scores, R = 200, target.with.covariates = FALSE)
qini <- maq(reward, cost, scores, R = 200)

# In some settings we may want to plot using one of R's many plot libraries.
# The plot method invisibly returns the plot data we can use for this purpose.
df.qini.baseline <- plot(qini.avg)
df.qini <- plot(qini, add = TRUE, col = 2)

# Make an alternate plot style, using, for example, ggplot.
ggplot(df.qini, aes(x = spend, y = gain)) +
  geom_ribbon(aes(ymin = gain - 1.96 * std.err,
                  ymax = gain + 1.96 * std.err),
              fill = "lightgray") +
  geom_line(linewidth = 2) +
  ylab("Policy value") +
  xlab("Fraction treated") +
  geom_line(data = df.qini.baseline, aes(x = spend, y = gain), lty = 2)
}
}

}
