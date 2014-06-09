#' Defunct functions/classes/methods in the MinimumDistance package
#'
#' The function, class, or data object you asked is defunct.
#'
#' @name coerce
#' @aliases coerce,RangedDataCNV,GRanges-class
#' @keywords internal
#' @rdname Defunct
#' @export
setAs("RangedDataCNV", "GRanges", function(from, to){
  .Defunct("RangedDataCNV class is defunct")
})
