## This is not general.  It assumes that we have (1) beadstudio type
## output and (2) that the output is in a format that can be handled
## by read.bsfiles.  This should be a function and not a method for
## class character.
setMethod("GenomeAnnotatedDataFrameFrom", signature(object="character"),
	  function(object, annotationPkg, genome, ...){
		  ##check if object is a file
		  if(!file.exists(object)) message("File ", object, " does not exist")
		  dat <- read.bsfiles(filenames=object)
		  GenomeAnnotatedDataFrameFrom(dat, annotationPkg=annotationPkg,
					       genome=genome, ...)
	  })

setMethod("sampleNames2", signature(object="AnnotatedDataFrame"),
	  function(object){
		  ## in order to allow duplicate fathers and mothers ids,
		  ## make.unique() was used to create the rownames for the annotated data frames.
		  originalNames(row.names(object@data))
	  })




