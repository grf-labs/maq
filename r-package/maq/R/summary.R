#' Qini curve summary.
#'
#' Get a data.frame with columns equal to \[B, Q(B), std.err(Q(B)), i, k\], where
#' i is the unit and k the treatment arm that is optimal to assign at a spend level B.
#'
#' @param object A maq object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A data.frame making up the elements of the estimated Qini curve.
#' @method summary maq
#' @export
summary.maq <- function(object,
                        ...) {

  data.frame(
    spend = object[["_path"]]$spend,
    gain = object[["_path"]]$gain,
    std.err = object[["_path"]]$std.err,
    unit.allocation = object[["_path"]]$ipath + 1, # +1: C++ index.
    arm.allocation = object[["_path"]]$kpath + 1
  )
}
