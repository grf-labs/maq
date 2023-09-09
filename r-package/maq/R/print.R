#' Print a maq object.
#' @param x A maq object.
#' @param ... Additional arguments (currently ignored).
#'
#' @method print maq
#' @export
print.maq <- function(x,
                      ...) {

  cat("MAQ object fit on", x$dim[1], "units and", x$dim[2], "arms with max budget", x$budget)
}
