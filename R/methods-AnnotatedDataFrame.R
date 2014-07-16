## This is not general.  It assumes that we have (1) beadstudio type
## output and (2) that the output is in a format that can be handled
## by read.bsfiles.  This should be a function and not a method for
## class character.

annotatedDataFrameFromArray <- function(object, byrow=FALSE, ...){
  if(dim(object)[[3]] > 0){
    object <- object[, , 1, drop=TRUE]
    ##res <- Biobase:::annotatedDataFrameFromMatrix(object, byrow=byrow, ...)
    res <- annotatedDataFrameFrom(object, byrow=byrow, ...)
  } else res <- annotatedDataFrameFrom(matrix(), byrow=byrow, ...)
  return(res)
}

setMethod("annotatedDataFrameFrom", signature(object="ff_array"),
          annotatedDataFrameFromArray)

setMethod("annotatedDataFrameFrom", signature(object="array"),
          annotatedDataFrameFromArray)

setMethod("GenomeAnnotatedDataFrameFrom", signature(object="character"),
	  function(object, annotationPkg, genome, ...){
		  ##check if object is a file
            if(!file.exists(object)) message("File ", object, " does not exist")
            dat <- read.bsfiles(object)
            id <- dat[[1]]
            i <- grep("Log.R", colnames(dat))
            dat <- GenomeAnnotatedDataFrameFrom(as.matrix(dat[[i]]),
                                                annotationPkg=annotationPkg,
                                                genome=genome, ...)
            rownames(dat) <- id
            dat
	  })

setMethod("sampleNames2", signature(object="AnnotatedDataFrame"),
	  function(object){
		  ## in order to allow duplicate fathers and mothers ids,
		  ## make.unique() was used to create the rownames for the annotated data frames.
            originalNames(row.names(object@data))
	  })
