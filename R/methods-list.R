setMethod("calculateMindist", signature(object="list"),
	  function(object, ...){
		  lapply(object, calculateMindist, ...)
	  })

setMethod("calculateMindist", signature(object="TrioSet"),
	  function(object, verbose=TRUE, ...){
        sns <- sampleNames(object)
	is.ff <- is(lrr(object), "ff")
	if(is.ff){
		invisible(open(lrr(object)))
	}
	md <- initializeBigMatrix("mindist", nr=nrow(object), nc=ncol(object),
				  vmode="double")
	if(verbose){
		message("\t\tComputing the minimum distance for ", ncol(object), " files.")
		pb <- txtProgressBar(min=0, max=ncol(object), style=3)
	}
	for(j in seq(length=ncol(object))){
		if(verbose) setTxtProgressBar(pb, j)
		LRR <- lrr(object)[, j, ]
		md[, j] <- calculateMindist(LRR)
	}
	if(verbose) close(pb)
	if(is.ff){
		close(md)
		close(lrr(object))
	}
	return(md)
})
