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

## this should probably be moved to VanillaICE, then imported by MinimumDistance
setMethod("GenomeAnnotatedDataFrameFrom", signature(object="array"),
	  function(object, annotationPkg){
		  GenomeAnnotatedDataFrameFromArray(object, annotationPkg)
	  })

GenomeAnnotatedDataFrameFromArray <- function(object, annotationPkg){
	## coerce to matrix
	dims <- dim(object)
	is.array <- length(dims) == 3
	if(is.array){
		res <- oligoClasses:::GenomeAnnotatedDataFrameFromMatrix(object[, , 1], annotationPkg)
	} else {
		##dim(object) <- dim(object)[c(1,2)]
		res <- oligoClasses:::GenomeAnnotatedDataFrameFromMatrix(object, annotationPkg)
	}
	res
}

## This is not general.  It assumes that we have (1) beadstudio type
## output and (2) that the output is in a format that can be handled
## by read.bsfiles.  This should be a function and not a method for
## class character.
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




