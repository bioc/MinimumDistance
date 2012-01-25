##---------------------------------------------------------------------------
##
## list
##
##---------------------------------------------------------------------------

setMethod("mad2", signature(object="list"),
	  function(object, byrow, pedigree, ...){
		  madList(object, byrow, pedigree, ...)
	  })

setMethod("mad2", signature(object="matrix"),
	  function(object, byrow, pedigree, ...){
		  madList(list(object), byrow, pedigree, ...)
	  })

setMethod("mad2", signature(object="array"),
	  function(object, byrow, pedigree, ...){
		  madList(list(object), byrow, pedigree, ...)
	  })

madList <- function(object, byrow, pedigree, ...){
	dims <- dim(object[[1]])
	if(length(dims) != 2 && length(dims) != 3)
		stop("Elements of list must be a matrix or an array")
	isff <- is(object[[1]], "ff")
	if(isff){
		lapply(object, open)
	}
	is.matrix <- ifelse(length(dims) == 2, TRUE, FALSE)
	if(!byrow){ ## by column
		if(is.matrix){
			mads <- madFromMatrixList(object, byrow=FALSE)
		} else { ## array
			## for parallelization, it would be better to
			## pass the ff object to the worker nodes,
			## calculate the mad, and return the mad.
			F <- lapply(object, function(x) as.matrix(x[, , 1]))
			M <- lapply(object, function(x) as.matrix(x[, , 2]))
			O <- lapply(object, function(x) as.matrix(x[, , 3]))
			mads.father <- madFromMatrixList(F, byrow=FALSE)
			mads.mother <- madFromMatrixList(M, byrow=FALSE)
			mads.offspr <- madFromMatrixList(O, byrow=FALSE)
			mads <- cbind(mads.father, mads.mother, mads.offspr)
			colnames(mads) <- c("F", "M", "O")
		}
	} else {## by row
		if(is.matrix){
			mads <- madFromMatrixList(object, byrow=FALSE)
		} else {
			if(ncol(object[[1]]) > 2){
			## for parallelization, it would be better to
			## pass the ff object to the worker nodes,
			## calculate the mad, and return the mad.
				stopifnot(!missing(pedigree))
				colindex <- which(!duplicated(fatherNames(pedigree)) & !duplicated(motherNames(pedigree)))
				##mindex <- which(!duplicated(motherNames(pedigree)))
				##F <- lapply(object, function(x) x[, findex, 1])
				##M <- lapply(object, function(x) x[, mindex, 2])
				O <- lapply(object, function(x) as.matrix(x[, colindex, 3]))
				##mads.f <- madFromMatrixList(F, byrow=TRUE)
				##mads.m <- madFromMatrixList(M, byrow=TRUE)
				mads <- madFromMatrixList(O, byrow=TRUE)
				names(mads) <- names(object)                                
				##mads <- cbind(mads.f, mads.m, mads.o)
				##colnames(mads) <- c("F", "M", "O")
			} else mads <- NULL
		}
	}
	if(isff) lapply(object, close)
	return(mads)
}


madFromMatrixList <- function(object, byrow=TRUE){
	if(!byrow){
		## this could be done more efficiently by following the
		## apply example in the foreach documentation...
		ilist <- splitIndicesByLength(seq_len(ncol(object[[1]])), 100)
		ispar <- !is.null(getCluster())
		if(ispar){
			Xlist <- foreach(i=ilist, .packages="MinimumDistance") %dopar% MinimumDistance:::stackListByColIndex(object, i)
			mads <- foreach(i = Xlist, .packages="MinimumDistance") %dopar% apply(i, 2, mad, na.rm=TRUE)
		} else {
			Xlist <- foreach(i=ilist, .packages="MinimumDistance") %do% MinimumDistance:::stackListByColIndex(object, i)
			mads <- foreach(i = Xlist, .packages="MinimumDistance") %do% apply(i, 2, mad, na.rm=TRUE)
		}
		mads <- unlist(mads)
		names(mads) <- colnames(object[[1]])
		mads
	} else {
		mads <- foreach(x = object, .packages="MinimumDistance") %do% VanillaICE:::rowMAD(x, na.rm=TRUE)
                if( !is.null(dim(mads[[1]])) & !is.null(rownames(object[[1]]))){
			labelrows <- function(x, fns) {
				rownames(x) <- fns
				return(x)
			}
			mads <- foreach(x=mads, fns=lapply(object, rownames)) %do% labelrows(x=x, fns=fns)
                }
	}
	return(mads)
}

##---------------------------------------------------------------------------
##
## TrioSet
##
##---------------------------------------------------------------------------

##setMethod("mad", signature(x="TrioSet"), function(x) x@mad)
##setMethod("mad.marker", signature(x="TrioSet"), function(x) {
##	res <- fData(x)$marker.mad
##	if(is.null(res)){
##		## message("Requesting NULL...too few samples to estimate sd")
##		##J <- ncol(x)
##		if(ncol(x) > 1) browser() ## which sample to use
##		tmp <- mad.sample(x)[, "O"]
##		res <- rep(tmp, nrow(x))
##	}
##	res
##})


##setMethod("mad.sample", signature(x="TrioSet"), function(x) x@mad)
##
##setReplaceMethod("mad.sample", signature(x="TrioSet", value="matrix"),
##	  function(x, value){
##		  x@mad <- value
##		  return(x)
##	  })
##
##setReplaceMethod("mad.marker", signature(x="TrioSet", value="numeric"),
##	  function(x, value){
##		  fData(x)$marker.mad <- value
##		  return(x)
##	  })
##
##setReplaceMethod("mad.marker", signature(x="TrioSet", value="matrix"),
##	  function(x, value){
##		  fData(x)$marker.mad <- value[, 3]
##		  fData(x)$marker.mad.father <- value[, 1]
##		  fData(x)$marker.mad.mother <- value[, 2]
##		  return(x)
##	  })
##
##setReplaceMethod("mad.marker", signature(x="TrioSet", value="NULL"),
##	  function(x, value){
##		  if(!is.null(value)){
##			  fData(x)$marker.mad <- value[, 3]
##			  fData(x)$marker.mad.father <- value[, 1]
##			  fData(x)$marker.mad.mother <- value[, 2]
##		  } else{
##			  x$marker.mad <- NULL
##		  }
##		  return(x)
##	  })
##
##setMethod("mad.mindist", signature(x="TrioSet"),
##	  function(x){
##		  x$mindist.mad
##	  })
##
##setReplaceMethod("mad.mindist", signature(x="TrioSet", value="numeric"),
##		 function(x, value){
##			 ## store in phenodata
##			 x$mindist.mad <- value
##			 return(x)
##		 })
##

##---------------------------------------------------------------------------
##
## TrioSetList
##
##---------------------------------------------------------------------------
##setMethod("mad", signature(x="TrioSetList"), function(x) mad(x[[1]]))

##setReplaceMethod("mad.sample", signature(x="TrioSetList", value="matrix"),
##		 function(x, value){
##			 for(i in seq_along(x)){
##				 mad.sample(x[[i]]) <- value
##			 }
##			 return(x)
##		 })

##setReplaceMethod("mad.mindist", signature(x="TrioSetList", value="numeric"),
##		 function(x, value){
##			 ## for genomewide estimation of mad parameter
##			 for(i in seq_along(x)){
##				 mad.mindist(x[[i]]) <- value
##			 }
##			 return(x)
##		 })
##
##setReplaceMethod("mad.mindist", signature(x="TrioSetList", value="list"),
##		 function(x, value){
##			 ## for chromosome-specific estimation of mad parameter
##			 for(i in seq_along(x)){
##				 mad.mindist(x[[i]]) <- value[[i]]
##			 }
##			 return(x)
##		 })
##
##setReplaceMethod("mad.marker", signature(x="TrioSetList", value="list"),
##		 function(x, value){
##			 for(i in seq_along(x)){
##				 mad.marker(x[[i]]) <- value[[i]]
##			 }
##			 return(x)
##		 })
##
##setReplaceMethod("mad.marker", signature(x="TrioSetList", value="NULL"),
##		 function(x, value){
##			 for(i in seq_along(x)){
##				 mad.marker(x[[i]]) <- NULL
##			 }
##			 return(x)
##		 })
##
##
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
