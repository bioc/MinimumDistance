#' Defunct functions/classes/methods in the MinimumDistance package
#'
#' The function, class, or data object you asked is defunct.
#'
#' @aliases coerce,RangedDataCNV,GRanges-class
#' @keywords internal
#' @rdname Defunct
#' @name Defunct
NULL
setAs("RangedDataCNV", "GRanges", function(from, to){
  .Defunct("RangedDataCNV class is defunct")
})
