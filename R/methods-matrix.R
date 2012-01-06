setMethod("calculateMindist", signature(object="array"),
	  function(object, ...){
		  ##stopifnot(ncol(object)==3)
		  d1 <- object[, , "O"] - object[, , "F"]
		  d2 <- object[, , "O"] - object[, , "M"]
		  I <- as.numeric(abs(d1) <= abs(d2))
		  md <- I*d1 + (1-I)*d2
		  md <- as.matrix(md)
		  colnames(md) <- colnames(object)
		  return(md)
	  })

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
