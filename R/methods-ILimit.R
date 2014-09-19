#' compute the range of an ILimit instance
#'
#' The range method for class ILimit is used internally in
#' MinimumDistance.
#' @param x a \code{ILimit} object
#' @param ... ignored
#' @param na.rm logical. If TRUE, missing values are removed.
setMethod("range", "ILimit", function(x, ..., na.rm=FALSE) c(start(x), end(x)))

setMethod("seq_along2", "ILimit", function(along.with){
  seq(start(along.with), end(along.with), 1)
})

ILimit <- function(...) as(IRanges(...), "ILimit")
