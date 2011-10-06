annotatedDataFrameFromArray <- function(object, byrow=FALSE, ...){
	if(dim(object)[[3]] > 0){
		object <- object[, , 1, drop=TRUE]
		res <- Biobase:::annotatedDataFrameFromMatrix(object, byrow=byrow, ...)
	} else res <- Biobase:::annotatedDataFrameFromMatrix(matrix(), byrow=byrow, ...)
	return(res)
}
setMethod("annotatedDataFrameFrom", signature(object="ff_array"),
	  annotatedDataFrameFromArray)
setMethod("annotatedDataFrameFrom", signature(object="array"),
	  annotatedDataFrameFromArray)
