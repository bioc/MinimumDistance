##---------------------------------------------------------------------------
##
## list
##
##---------------------------------------------------------------------------
setMethod("mad2", signature(object="list"),
	  function(object, byrow, ...){
		  is.matrix <- is(object[[1]], "matrix") || is(object[[1]], "ff_matrix")
		  is.array <- class(object[[1]])=="array" || is(object[[1]], "ff_array")
		  ##is.matrix <- is(object[[1]], "matrix")
		  ##is.array <- is(object[[1]], "array")
		  stopifnot(is.matrix || is.array)
		  if(!byrow){ ## by column
			  if(is.matrix){
				  mads <- rep(NA, ncol(object[[1]]))
				  names(mads) <- colnames(object[[1]])
				  ## to avoid memory problems
				  ilist <- splitIndicesByLength(seq_len(ncol(object[[1]])), 100)
				  for(i in seq_along(ilist)){
					  index <- ilist[[i]]
					  X <- stackListByColIndex(object, index)
					  mads[index] <- apply(X, 2, mad, na.rm=TRUE)
				  }
			  } else {
				  J <- dim(object[[1]])[[3]]
				  mads <- matrix(NA, ncol(object[[1]]), 3)
				  colnames(mads) <- c("F", "M", "O")
				  rownames(mads) <- colnames(object[[1]])
				  ## to avoid memory problems
				  ilist <- splitIndicesByLength(seq_len(ncol(object[[1]])), 100)
				  for(i in seq_along(ilist)){
					  index <- ilist[[i]]
					  for(j in seq_len(J)){
						  X <- stackListByColIndex(object, index, j)
						  mads[index, j] <- apply(X, 2, mad, na.rm=TRUE)
					  }
				  }
			  }
		  } else {
			  mads <- vector("list", length(object))
			  if(is.matrix){
				  for(i in seq_along(object)){
					  ##j <- indexlist[[i]]
					  mads[[i]] <- rowMAD(object[[i]], na.rm=TRUE)
					  rownames(mads[[i]]) <- rownames(object[[i]])
				  }
			  } else {
				  J <- dim(object[[1]])[[3]]
				  ##mads <- matrix(NA, nr, 3)
				  ##dimnames(mads) <- list(fns, c("F", "M", "O"))
				  for(i in seq_along(object)){
					  mads[[i]] <- matrix(NA, nrow(object[[i]]), 3)
					  dimnames(mads[[i]]) <- list(rownames(object[[i]]),
								      c("F","M", "O"))
					  for(j in seq_len(J)){
						  mads[[i]][, j] <- rowMAD(object[[i]][, ,j], na.rm=TRUE)
					  }
				  }
			  }
		  }
		  return(mads)
	  })

##---------------------------------------------------------------------------
##
## TrioSet
##
##---------------------------------------------------------------------------

setMethod("mad", signature(x="TrioSet"), function(x) x@mad)
setMethod("mad.marker", signature(x="TrioSet"), function(x) fData(x)$marker.mad)
setMethod("mad.sample", signature(x="TrioSet"), function(x) x@mad)

setReplaceMethod("mad.sample", signature(x="TrioSet", value="matrix"),
	  function(x, value){
		  x@mad <- value
		  return(x)
	  })

setReplaceMethod("mad.marker", signature(x="TrioSet", value="numeric"),
	  function(x, value){
		  fData(x)$marker.mad <- value
		  return(x)
	  })

setReplaceMethod("mad.marker", signature(x="TrioSet", value="matrix"),
	  function(x, value){
		  fData(x)$marker.mad <- value[, 3]
		  fData(x)$marker.mad.father <- value[, 1]
		  fData(x)$marker.mad.mother <- value[, 2]
		  return(x)
	  })

setMethod("mad.mindist", signature(x="TrioSet"),
	  function(x){
		  x$mindist.mad
	  })

setReplaceMethod("mad.mindist", signature(x="TrioSet", value="numeric"),
		 function(x, value){
			 ## store in phenodata
			 x$mindist.mad <- value
			 return(x)
		 })


##---------------------------------------------------------------------------
##
## TrioSetList
##
##---------------------------------------------------------------------------
setMethod("mad", signature(x="TrioSetList"), function(x) mad(x[[1]]))

setReplaceMethod("mad.sample", signature(x="TrioSetList", value="matrix"),
		 function(x, value){
			 for(i in seq_along(x)){
				 mad.sample(x[[i]]) <- value
			 }
			 return(x)
		 })

setReplaceMethod("mad.mindist", signature(x="TrioSetList", value="numeric"),
		 function(x, value){
			 ## for genomewide estimation of mad parameter
			 for(i in seq_along(x)){
				 mad.mindist(x[[i]]) <- value
			 }
			 return(x)
		 })

setReplaceMethod("mad.mindist", signature(x="TrioSetList", value="list"),
		 function(x, value){
			 ## for chromosome-specific estimation of mad parameter
			 for(i in seq_along(x)){
				 mad.mindist(x[[i]]) <- value[[i]]
			 }
			 return(x)
		 })

setReplaceMethod("mad.marker", signature(x="TrioSetList", value="list"),
		 function(x, value){
			 for(i in seq_along(x)){
				 mad.marker(x[[i]]) <- value[[i]]
			 }
			 return(x)
		 })

##setReplaceMethod("mad.marker", signature(x="TrioSetList", value="matrix"),
##		 function(x, value){
##			 nr <- sum(sapply(x, nrow))
##			 indexlist <- split(seq_len(nr), rep(seq_along(x), sapply(x, nrow)))
##			 for(i in seq_along(x)){
##				 j <- indexlist[[i]]
##				 mad.marker(x[[i]]) <- value[j, ]
##			 }
##			 return(x)
##		 })
