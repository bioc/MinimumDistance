setMethod("calculateMindist", signature(object="arrayORff_array"),
	  function(object, outdir, ...){
		  ##stopifnot(ncol(object)==3)
		  calculateMindistFromArray(object, outdir, ...)
	  })

calculateMindistFromArray <- function(object, outdir=ldPath(), ...){
	isff <- is(object, "ff")
	if(is.null(getCluster())) registerDoSEQ()
	if(isff){
		require("ff")
		## so that the worker nodes put the ff objects in the same directory
		ldPath(outdir)
		md <- initializeBigMatrix("mindist", nr=nrow(object), nc=ncol(object), vmode="double")
		for(j in seq_len(ncol(object))){
			d1 <- object[, j, 3] - object[, j, 1] ## offspring - father
			d2 <- object[, j, 3] - object[, j, 2] ## offspring - mother
			I <- as.numeric(abs(d1) <= abs(d2))
			md[, j] <- I*d1 + (1-I)*d2
		}
		colnames(md) <- colnames(object)
	} else {
		d1 <- object[, , 3] - object[, , 1] ##offspring - father
		d2 <- object[, , 3] - object[, , 2] ##offspring - mother
		I <- as.numeric(abs(d1) <= abs(d2))
		md <- I*d1 + (1-I)*d2
		md <- as.matrix(md)
		colnames(md) <- colnames(object)
	}
	return(md)
}

##setMethod("cnEmission", signature(object="array"),
##	  function(object, stdev, k=5, cnStates, is.log, is.snp,
##		   normalIndex, verbose=TRUE, ...){
##		  emit <- array(NA, dim=c(dim(object), length(cnStates)))
##		  for(j in seq_len(ncol(object))){
##			  emit[, j, , ] <- cnEmission(object[, j, ],
##						    stdev,
##						    k,
##						    cnStates,
##						    is.log,
##						    is.snp,
##						    normalIndex,
##						    verbose, ...)
##		  }
##		  return(emit)
##	  })
##
##setMethod("bafEmission", signature(object="array"),
##	  function(object, is.snp, prOutlier=1e-3, p.hom=0.95, ...){
##		  emit <- array(NA, dim=c(dim(object), 6))
 ##		  for(j in seq_len(ncol(object))){
##			  emit[, j, , ] <- bafEmission(object[, j, ],
##						       is.snp,
##						       prOutlier,
##						       p.hom)
##		  }
##		  return(emit)
##	  })
