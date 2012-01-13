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

setMethod("GenomeAnnotatedDataFrameFrom", signature(object="array"),
	  function(object, annotationPkg){
		  GenomeAnnotatedDataFrameFromArray(object, annotationPkg)
	  })

GenomeAnnotatedDataFrameFromArray <- function(object, annotationPkg){
	## coerce to matrix
	if(nrow(object) > 1 & ncol(object) > 1){
		oligoClasses:::GenomeAnnotatedDataFrameFromMatrix(object[, , 1], annotationPkg)
	} else {
		dim(object) <- dim(object)[c(1,2)]
		oligoClasses:::GenomeAnnotatedDataFrameFromMatrix(object, annotationPkg)
	}
}

setMethod("GenomeAnnotatedDataFrameFrom", signature(object="character"),
	  function(object, annotationPkg){
		  ##check if object is a file
		  if(!file.exists(object)) message("File ", object, " does not exist")
		  dat <- read.bsfiles(filenames=object)
		  GenomeAnnotatedDataFrameFrom(dat, annotationPkg)
	  })

setMethod("sampleNames2", signature(object="AnnotatedDataFrame"),
	  function(object){
		  ## in order to allow duplicate fathers and mothers ids,
		  ## make.unique() was used to create the rownames for the annotated data frames.
		  originalNames(row.names(object@data))
	  })




