##---------------------------------------------------------------------------
##
## list
##
##---------------------------------------------------------------------------

#' @param pedigree an object of class \code{Pedigree}
#' @aliases mad2,list-method
#' @rdname mad2
setMethod("mad2", signature(object="list"),
	  function(object, byrow, pedigree, ...){
		  madList(object, byrow, pedigree, ...)
	  })

#' @aliases mad2,TrioSetList-method
#' @rdname mad2
setMethod("mad2", signature(object="TrioSetList"),
	  function(object, byrow, ...){
		  madTrioSetList(object, byrow)
	  })

madTrioSetList <- function(object, byrow){
	madList(lrr(object), byrow=byrow, pedigree=pedigree(object))
}

#' @aliases mad2,matrix-method
#' @rdname mad2
setMethod("mad2", signature(object="matrix"),
	  function(object, byrow, pedigree, ...){
		  madList(list(object), byrow, pedigree, ...)
	  })

#' @aliases mad2,array-method
#' @rdname mad2
setMethod("mad2", signature(object="array"),
	  function(object, byrow, pedigree, ...){
		  madList(list(object), byrow, pedigree, ...)
	  })

madList <- function(object, byrow, pedigree, ...){
	dims <- dim(object[[1]])
	if(length(dims) != 2 && length(dims) != 3)
		stop("Elements of list must be a matrix or an array")
	isff <- is(object[[1]], "ff")
	if(isff) lapply(object, open)
	is.matrix <- ifelse(length(dims) == 2, TRUE, FALSE)
	if(!byrow){ ## by column
		if(is.matrix){
			mads <- madFromMatrixList(object, byrow=FALSE)
		} else { ## array
			## for parallelization, it would be better to
			## pass the ff object to the worker nodes,
			## calculate the mad, and return the mad.
			if(!isff){
				F <- lapply(object, function(x) as.matrix(x[, , 1]))
				M <- lapply(object, function(x) as.matrix(x[, , 2]))
				O <- lapply(object, function(x) as.matrix(x[, , 3]))
				mads.father <- madFromMatrixList(F, byrow=FALSE)
				mads.mother <- madFromMatrixList(M, byrow=FALSE)
				mads.offspr <- madFromMatrixList(O, byrow=FALSE)
				if(!missing(pedigree)){
					names(mads.father) <- fatherNames(pedigree)
					names(mads.mother) <- motherNames(pedigree)
					names(mads.offspr) <- offspringNames(pedigree)
					mads <- data.frame(F=I(mads.father),
							   M=I(mads.mother),
							   O=I(mads.offspr))
				} else {
					mads <- cbind(mads.father, mads.mother, mads.offspr)
					colnames(mads) <- c("F", "M", "O")
				}
			} else { ## big data
				madForFFmatrix <- function(xlist, i, j){
					## j=1 father
					## j=2 mother
					## j=3 offspring
					res <- lapply(xlist, function(x, i, j){
						m <- as.matrix(x[, i, j])
					}, i=i, j=j)
					res <- do.call("rbind", res)/100
					mads <- apply(res, 2, mad, na.rm=TRUE)
					mads
				}
				indexList <- splitIndicesByLength(seq_len(ncol(object[[1]])), 50)
				i <- NULL
				mads.father <- foreach(i=indexList, .combine="c") %do% madForFFmatrix(xlist=object, i=i, j=1)
				mads.mother <- foreach(i=indexList, .combine="c") %do% madForFFmatrix(xlist=object, i=i, j=2)
				mads.offspr <- foreach(i=indexList, .combine="c") %do% madForFFmatrix(xlist=object, i=i, j=3)
				if(!missing(pedigree)){
					names(mads.father) <- fatherNames(pedigree)
					names(mads.mother) <- motherNames(pedigree)
					names(mads.offspr) <- offspringNames(pedigree)
					mads <- data.frame(F=I(mads.father),
							   M=I(mads.mother),
							   O=I(mads.offspr))
				} else {
					mads <- cbind(mads.father, mads.mother, mads.offspr)
					colnames(mads) <- c("F", "M", "O")
				}
			}
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
			} else {
				warning("Too few samples to calculate across sample variance. Returning NULL.")
				mads <- NULL
			}
		}
	}
	if(isff) lapply(object, close)
	return(mads)
}

madFromMatrixList <- function(object, byrow=TRUE){
  if(isPackageLoaded("ff")) pkgs <- c("ff", "MinimumDistance") else pkgs <- "MinimumDistance"
  if(!byrow){
    ## this could be done more efficiently by following the
    ## apply example in the foreach documentation...
    ilist <- splitIndicesByLength(seq_len(ncol(object[[1]])), 100)
    i <- NULL
    Xlist <- foreach(i=ilist, .packages=pkgs) %dopar% stackListByColIndex(object, i)
    mads <- foreach(i = Xlist, .packages=pkgs) %dopar% apply(i/100, 2, mad, na.rm=TRUE)
    mads <- unlist(mads)
    names(mads) <- colnames(object[[1]])
    mads
  } else {
    x <- NULL
    mads <- foreach(x = object, .packages=pkgs) %do% rowMAD(x/100, na.rm=TRUE)
    if( !is.null(dim(mads[[1]])) & !is.null(rownames(object[[1]]))){
      labelrows <- function(x, fns) {
        rownames(x) <- fns
        return(x)
      }
      fns <- NULL
      mads <- foreach(x=mads, fns=lapply(object, rownames)) %do% labelrows(x=x, fns=fns)
    }
  }
  return(mads)
}
