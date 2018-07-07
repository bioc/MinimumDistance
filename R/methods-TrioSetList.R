setMethod("initialize", signature(.Object="TrioSetList"),
	  function(.Object,
		   pedigreeData=Pedigree(),
		   assayDataList=AssayDataList(BAF=BAF, logRRatio=logRRatio),
		   logRRatio=array(NA, dim=c(0,0,3)),
		   BAF=array(NA, dim=dim(logRRatio)),
		   featureDataList=GenomeAnnotatedDataFrameFromList(assayDataList),
		   chromosome=integer(),
		   phenoData=annotatedDataFrameFrom(assayDataList, byrow=FALSE),
		   fatherPhenoData=annotatedDataFrameFrom(assayDataList, byrow=FALSE),
		   motherPhenoData=annotatedDataFrameFrom(assayDataList, byrow=FALSE),
		   genome=c("hg19", "hg18"),
		   ...){
		  .Object <- callNextMethod(.Object,
					    pedigree=pedigreeData,
					    assayDataList=assayDataList,
					    featureDataList=featureDataList,
					    phenoData=phenoData,
					    fatherPhenoData=fatherPhenoData,
					    motherPhenoData=motherPhenoData,
					    chromosome=chromosome,
					    genome=match.arg(genome),
					    ...)
		  if(all(sapply(featureDataList, nrow) == 0)) .Object@genome <- ""
		  .Object
	  })

setMethod("elementNROWS", signature(x="TrioSetList"), function(x){
	if(length(x) == 0) return(0L)
	as.integer(sapply(featureData(x), nrow))
})

setValidity("TrioSetList", function(object){
	l <- elementNROWS(object)
	if(any(l > 0)){
		if(!genomeBuild(object) %in% c("hg19", "hg18"))
			return(paste("genome is ", genomeBuild(object), ", but must be 'hg18' or 'hg19'.", sep=""))
	}
	nms <- ls(assayData(object))
	if(!all(c("BAF", "logRRatio") %in% nms)){
		msg <- "BAF and logRRatio are required elements of the assayData"
		return(msg)
	}
	if(length(object) > 0){
		msg <- validAssayDataDims(assayData(object))
		if(!all(msg == TRUE)) return(msg)
		elt <- (ls(assayDataList(object)))[[1]]
		b <- assayDataList(object)[[elt]]
		if(length(chromosome(object)) != length(b)){
			return("chromosome slot must be the same length as the length of the list for each assayData element")
		}
	}
	validObject(pedigree(object))
	if(!identical(sampleNames(object), sampleNames(phenoData(object)))){
		stop("sampleNames of TrioSetList object must be the same as the sampleNames of the phenoData")
	}
	if(!identical(fatherNames(object), originalNames(sampleNames(fatherPhenoData(object))))){
		stop("fatherNames of TrioSetList object must be the same as the sampleNames of the fatherPhenoData")
	}
	if(!identical(motherNames(object), originalNames(sampleNames(motherPhenoData(object))))){
		stop("motherNames of TrioSetList object must be the same as the sampleNames of the motherPhenoData")
	}
	if(length(featureDataList(object)) != length(chromosome(object))){
		return("each chromosome should have an element in the featureDataList")
	}
	if(length(featureDataList(object)) > 0){
		isGAD <- sapply(featureDataList(object), function(x) is(x, "GenomeAnnotatedDataFrame"))
		if(!all(isGAD)) return("featureDataList must be comprised of GenomeAnnotatedDataFrame(s)")
	}
})

setMethod("updateObject", signature(object="TrioSetList"),
	  function(object, ..., verbose=FALSE){
		  if (verbose) message("updateObject(object = 'TrioSetList')")
		  if(!is(object@featureDataList[[1]], "GenomeAnnotatedDataFrame")){
			  fdlist <- lapply(object@featureDataList, updateObject)
			  object@featureDataList <- fdlist
		  }
		  return(object)
	  })

##setMethod("lapply", signature(X="TrioSetList"),
##	  function(X, FUN, ...){
##		  res <- vector("list", length(X))
##		  for(i in seq_along(X)){
##			  res[[i]] <- FUN(X[[i]], ...)
##		  }
##		  res <- new("TrioSetList",
##			     assayData=
##
##		  return(res)
##	  })


GenomeAnnotatedDataFrameFromList <- function(object, annotationPkg){
	nms <- ls(object)
	elt <- object[[nms[1]]]
	fdlist <- vector("list", length(elt))
	for(i in seq_along(elt)){
		##fdlist[[i]] <- GenomeAnnotatedDataFrameFromArray(elt[[i]], annotationPkg)
		fdlist[[i]] <- GenomeAnnotatedDataFrameFrom(elt[[i]], annotationPkg)
	}
	return(fdlist)
}


#' Constructor for \code{TrioSetList} class
#'
#' The \code{TrioSetList} class has been deprecated and may be removed in
#' a future release. Use \code{MinDistExperiment} instead.
#'
#' @param chromosome integer vector of chromosome names
#' @param pedigreeData a \code{Pedigree} object
#' @param sample.sheet a \code{data.frame} containing sample covariates
#' @param row.names a character vector
#' @param lrr a matrix of log R ratios
#' @param baf a matrix of B allele frequencies
#' @param featureData a \code{GenomeAnnotatedDataFrame}
#' @param cdfname a character string indicating the annotation package
#' @param ffname prefix for ff-filenames
#' @param genome character string indicating genome build
#' @export
TrioSetList <- function(chromosome=integer(),
			pedigreeData=Pedigree(),
			sample.sheet,
			row.names=NULL,
			lrr, baf,
			featureData,
			cdfname,
			ffname="",
			genome){
  if(!missing(lrr)){
    if(!is(lrr[1,1], "integer")){
      stop("lrr should be a matrix of integers. Use integerMatrix(x, scale=100) for the transformation")
    }
    if(!is(baf[1,1], "integer")){
      stop("baf should be a matrix of integers.  Use integerMatrix(x, scale=1000) for the transformation")
    }
    if(missing(genome)) stop("Argument genome is missing.  Must specify UCSC genome build genome ('hg18' or 'hg19').")
    if(!missing(pedigreeData)){
      if(!all(fatherNames(pedigreeData) %in% colnames(lrr))) stop("column names of lrr and baf matrices must match names the pedigree file")
      if(!all(motherNames(pedigreeData) %in% colnames(lrr))) stop("column names of lrr and baf matrices must match names the pedigree file")
      if(!all(offspringNames(pedigreeData) %in% colnames(lrr))) stop("column names of lrr and baf matrices must match names the pedigree file")
    }
  }
  if(nrow(pedigreeData) > 0 & !(missing(lrr) | missing(baf))){
    if(!missing(sample.sheet)){
      if(is.null(row.names)){
        row.names <- rownames(sample.sheet)
      }
      index <- row.names %in% allNames(pedigreeData)
      sample.sheet <- sample.sheet[index, ]
      row.names <- row.names[index]
      if(!all(row.names %in% allNames(pedigreeData))){
        stop("There are row.names for sample.sheet not in the pedigree object")
      }
      phenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE,
                                          sample.sheet=sample.sheet,
                                          which="offspring",
                                          row.names=row.names)
      fatherPhenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE,
                                                sample.sheet=sample.sheet,
                                                which="father",
                                                row.names=row.names)
      motherPhenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE,
                                                sample.sheet=sample.sheet,
                                                which="mother",
                                                row.names=row.names)
    }  else {
      phenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE, which="offspring")
      fatherPhenoData <- annotatedDataFrameFrom(pedigreeData, FALSE, which="father")
      motherPhenoData <- annotatedDataFrameFrom(pedigreeData, FALSE, which="mother")
    }
  }
  if(length(chromosome) > 0){
    if(!all(chromosome %in% 1:22)){
      stop("Only autosomal chromosomes (1, 2, ... , 22) allowed")
    }
    if(any(duplicated(chromosome)))
      stop("duplicated chromosomes present")
  } else {
    if(missing(lrr) & missing(baf))
      return(new("TrioSetList"))
  }
  if(missing(lrr) | missing(baf)){
    lrrlist <- baflist <- lapply(chromosome, function(x) array(NA, dim=c(0,0,3)))
    ad <- AssayDataList(BAF=baflist, logRRatio=lrrlist)
    object <- new("TrioSetList",
                  assayDataList=ad,
                  chromosome=chromosome)
    return(object)
  }
  if(!identical(rownames(lrr), rownames(baf))) stop("rownames of lrr and baf must be identical")
  if(missing(featureData)){
    if(missing(cdfname)) stop("if featureData is not supplied, a valid cdfname must be provided for annotating the markers")
    if(any(is.na(rownames(lrr)))){
      message("Removing rows with NA identifiers from lrr & baf matrices")
      lrr <- lrr[!is.na(rownames(lrr)), ]
      baf <- baf[!is.na(rownames(baf)), ]
    }
    ##featureData <- oligoClasses:::featureDataFrom(cdfname)
    featureData <- GenomeAnnotatedDataFrameFrom(lrr, cdfname, genome=genome)
    fD <- featureData[order(chromosome(featureData), position(featureData)), ]
    rm(featureData); gc()
  } else {
    if(!is(featureData, "GenomeAnnotatedDataFrame")) stop("featureData must be a GenomeAnnotatedDataFrame")
    fD <- featureData
  }
  if(length(chromosome) > 0){
    fD <- fD[fD$chromosome%in%chromosome, ]
  }
  if(!is.null(rownames(lrr))){
    is.present <- featureNames(fD) %in% rownames(lrr)
    if(!all(is.present)) fD <- fD[is.present, ]
    index <- match(featureNames(fD), rownames(lrr))
    lrr <- lrr[index, ]
    baf <- baf[index, ]
    if(!all(identical(rownames(lrr), sampleNames(fD))))
      stop("rownames of lrr must be the same as the featureNames for the featureData")
  }
  marker.list <- split(seq_along(sampleNames(fD)), fD$chromosome)
  np <- nrow(trios(pedigreeData))
  trio.names <- array(NA, dim=c(length(offspringNames(pedigreeData)), 1, 3))
  dimnames(trio.names) <- list(offspringNames(pedigreeData), "sampleNames", c("F", "M", "O"))
  trio.names[, "sampleNames", ] <- as.matrix(trios(pedigreeData))
  father.names <- originalNames(fatherNames(pedigreeData))
  mother.names <- originalNames(motherNames(pedigreeData))
  offspring.names <- offspringNames(pedigreeData)
  father.index <- match(father.names, colnames(lrr))
  mother.index <- match(mother.names, colnames(lrr))
  offspring.index <- match(offspring.names, colnames(lrr))
  chromosome <- unique(chromosome(fD))
  fdlist <- baflist <- lrrlist <- vector("list", length(chromosome))
  if(isPackageLoaded("ff")){
    if(ffname!=""){
      bafname <- paste(ffname, "baf", sep="_")
    } else bafname <- "baf"
    if(ffname!=""){
      lrrname <- paste(ffname, "lrr", sep="_")
    } else lrrname <- "lrr"
  }
  dns <- list(sampleNames(pedigreeData), c("F", "M", "O"))
  for(i in seq_along(marker.list)){
    ## Use the name of the offspring as the name for the trio:
    j <- marker.list[[i]]
    nr <- length(j)
    bafArray <- initializeBigArray(bafname, dim=c(nr, np, 3), vmode="integer")
    logRArray <- initializeBigArray(lrrname, dim=c(nr, np, 3), vmode="integer")
    dimnames(logRArray)[c(2,3)] <- dimnames(bafArray)[c(2,3)] <- dns
    logRArray[,,"F"] <- lrr[j, father.index]
    logRArray[,,"M"] <- lrr[j, mother.index]
    logRArray[,,"O"] <- lrr[j, offspring.index]
    bafArray[,,"F"] <- baf[j, father.index]
    bafArray[,,"M"] <- baf[j, mother.index]
    bafArray[,,"O"] <- baf[j, offspring.index]
    ## For each chromosome, create a TrioSet
    lrrlist[[i]] <- logRArray
    baflist[[i]] <- bafArray
    fdlist[[i]] <- fD[j, ]
  }
  ad <- AssayDataList(logRRatio=lrrlist,
                      BAF=baflist)
  object <- new("TrioSetList",
                assayDataList=ad,
                featureDataList=fdlist,
                chromosome=chromosome,
                pedigree=pedigreeData,
                fatherPhenoData=fatherPhenoData,
                motherPhenoData=motherPhenoData,
                phenoData=phenoData,
                genome=genome)
  return(object)
}



setMethod("featureNames", signature(object="TrioSetList"),
	  function(object){
		  lapply(featureDataList(object), sampleNames)
	  })

setMethod("position", signature(object="TrioSetList"),
	  function(object){
		  lapply(featureDataList(object), position)
	  })

setMethod("isSnp", signature(object="TrioSetList"),
	  function(object){
		  lapply(featureDataList(object), function(x) isSnp)
	  })

setMethod("allNames", signature(object="TrioSetList"), function(object) allNames(pedigree(object)))

#' @param object a \code{TrioSetList} object
#' @aliases pedigree,TrioSetList-method
#' @rdname TrioSetList-class
setMethod("pedigree", signature(object="TrioSetList"), function(object) object@pedigree)

#' @aliases trios,TrioSetList-method
#' @rdname TrioSetList-class
setMethod("trios", signature(object="TrioSetList"), function(object) trios(pedigree(object)))
setMethod("sampleNames", signature(object="TrioSetList"),
	  function(object) sampleNames(pedigree(object)))
setMethod("nrow", signature(x="TrioSetList"),
	  function(x){
	  sum(sapply(x, nrow))
  })
setMethod("ncol", signature(x="TrioSetList"),
	  function(x) ncol(x[[1]]))

#' @aliases offspringNames,TrioSetList-method
#' @rdname TrioSetList-class
setMethod("offspringNames", signature(object="TrioSetList"), function(object){
  offspringNames(pedigree(object))
})

setMethod("fatherNames", signature(object="TrioSetList"), function(object){
	fatherNames(pedigree(object))
})
setMethod("motherNames", signature(object="TrioSetList"), function(object){
	motherNames(pedigree(object))
})

setMethod("annotation", signature(object="TrioSetList"), function(object){
	annotation(object[[1]])
})

setMethod("dims", signature(x="TrioSetList"), function(x){
	nr <- nrow(x)
	nchr <- length(chromosome(x))
	ntrios <- ncol(baf(x)[[1]])
	dm <- c(nchr, ntrios, nr)
	names(dm) <- c("chromosomes", "trios", "features")
	return(dm)
})




setMethod("sampleNames", signature(object="TrioSetList"),
	  function(object) offspringNames(object))
##setReplaceMethod("sampleNames", signature(object="TrioSetList", value="character"),
##		 function(object, value){
##			 object <- lapply(object, function(x, value ){
##				 sampleNames(x) <- value
##				 return(x)
##				 }, value=value)
##			 object <- as(object, "TrioSetList")
##			 return(object)
##	 })

##setReplaceMethod(mindist, c("TrioSetList", "list"),
##		 function(object, value){
##
##})

setMethod("prune", signature(object="TrioSetList", ranges="RangedDataCNV"),
	  function(object, ranges, id, lambda, min.change, min.coverage,
		   scale.exp, verbose, ...){
		  rdList <- lapply(object, prune, ranges=ranges,
				   id=id,
				   lambda=lambda,
				   min.change=min.change,
				   min.coverage=min.coverage,
				   scale.exp=scale.exp,
				   verbose=verbose, ...)
		  return(rdList)
	  })


setMethod("assayData", signature(object="TrioSetList"),
	  function(object) assayDataList(object))
setMethod("storageMode", "TrioSetList", function(object) storageMode(assayData(object)))

setMethod("phenoData", signature(object="TrioSetList"),
	  function(object) object@phenoData)
setMethod("offspringPhenoData", signature(object="TrioSetList"),
	  function(object) phenoData(object))
setMethod("fatherPhenoData", signature(object="TrioSetList"),
	  function(object) object@fatherPhenoData)
setMethod("motherPhenoData", signature(object="TrioSetList"),
	  function(object) object@motherPhenoData)

setReplaceMethod("assayData", signature=signature(object="TrioSetList",
			      value="AssayData"),
                 function(object, value) {
			 object@assayDataList <- value
			 object
                 })

setReplaceMethod("phenoData", signature=signature(object="TrioSetList",
			      value="AnnotatedDataFrame"),
                 function(object, value) {
			 object@phenoData <- value
			 object
                       })

#' @param x a \code{TrioSetList}
#' @param i a numeric vector for subsetting the chromosomes  (optional)
#' @param j a numeric vector for subsetting trios (optional)
#' @param ... additional arguments passed to subsetting methods for matrices and data frames
#' @param drop logical. Whether to simplify matrices to numeric
#' vectors.  This should be left as FALSE.
#' @aliases "[",TrioSetList,ANY-method
#' @rdname TrioSetList-class
setMethod("[", signature(x="TrioSetList"),
	  function(x, i, j, ..., drop=FALSE){
		  ## using 'i' to subset markers does not really make
		  ## sense
		  ##
		  ## Use i to subset the list. example, x[1] is still a TrioSetList, but is one chromosome
		  ##
		  if(!missing(i) & !missing(j)){
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- ad[[elt]][i]
				  tmp <- lapply(tmp, function(x, j) {
					  x[, j, , drop=FALSE]
				  }, j=j)
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@chromosome <- chromosome(x)[i]
			  x@featureDataList <- featureDataList(x)[i]
			  x@pedigree <- pedigree(x)[j, ]
			  x@phenoData <- phenoData(x)[j, ]
			  x@fatherPhenoData <- fatherPhenoData(x)[j, ]
			  x@motherPhenoData <- motherPhenoData(x)[j, ]
		  }
		  if(!missing(i) & missing(j)){
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- ad[[elt]][i]
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@chromosome <- chromosome(x)[i]
			  x@featureDataList <- featureDataList(x)[i]
		  }
		  if(missing(i) & !missing(j)){
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- lapply(ad[[elt]], function(x, j) {
					  x[, j, , drop=FALSE]
				  }, j=j)
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@pedigree <- pedigree(x)[j, ]
			  x@phenoData <- phenoData(x)[j, ]
			  x@fatherPhenoData <- fatherPhenoData(x)[j, ]
			  x@motherPhenoData <- motherPhenoData(x)[j, ]
		  }
		  return(x)
	  })

#' @param exact ignored
#' @aliases "[[",TrioSetList,ANY,ANY-method
#' @rdname TrioSetList-class
setMethod("[[", signature(x="TrioSetList"),
	  function(x, i, j, ..., exact=TRUE){
		  if(missing(i)) return(x)
		  if(length(i) == 1){
			  lrrs <- lrr(x)[[i]]
			  bafs <- baf(x)[[i]]
			  fdlist <- featureDataList(x)[[i]]
			  x <- new("TrioSet",
				   logRRatio=lrrs,
				   BAF=bafs,
				   phenoData=phenoData(x),
				   fatherPhenoData=fatherPhenoData(x),
				   motherPhenoData=motherPhenoData(x),
				   pedigree=pedigree(x),
				   featureData=featureDataList(x)[[i]],
				   genome=genomeBuild(x))
		  } else {
			  stop("subscript out of bounds")
		  }
	  })

#' @aliases show,TrioSetList-method
#' @rdname TrioSetList-class
setMethod("show", signature(object="TrioSetList"),
	  function(object){
		  lo <- length(lrr(object))
		  cat(class(object), " of length ", lo, "\n", sep="")
		  cat("genome:", genomeBuild(object), "\n")
	  })

#' @aliases length,TrioSetList-method
#' @rdname TrioSetList-class
setMethod("length", signature(x="TrioSetList"), function(x) length(x@chromosome))


#' @aliases calculateMindist,TrioSetList-method
#' @rdname calculateMindist
setMethod("calculateMindist", signature(object="TrioSetList"),
	  function(object){
		  AssayDataList(calculateMindist(lrr(object)))
                })



setMethod("assayDataList", signature(object="TrioSetList"),
	  function(object)  object@assayDataList)

setMethod("featureDataList", signature(object="TrioSetList"),
	  function(object)  object@featureDataList)

setMethod("featureData", signature(object="TrioSetList"),
	  function(object)  object@featureDataList)

setMethod("lrr", signature(object="TrioSetList"),
	  function(object){
		  ##lapply(object, lrr)
		  assayDataList(object)[["logRRatio"]]
	  })

setMethod("baf", signature(object="TrioSetList"),
	  function(object){
		  ##lapply(object, baf)
		  assayDataList(object)[["BAF"]]
	  })

setMethod("chromosome", signature(object="TrioSetList"),
	  function(object, as.list=FALSE, ...){
		  ##lapply(object, chromosome)
		  if(!as.list) object@chromosome else chromosomeList(object)
	  })

setMethod("chromosomeList", signature(object="TrioSetList"),
	  function(object){
		  ##lapply(object, chromosome)
		  lrrs <- lrr(object)
		  chrom <- rep(object@chromosome, sapply(lrrs, nrow))
		  split(chrom, chrom)
	  })

setMethod("checkOrder", signature(object="TrioSetList"),
	  function(object, verbose=FALSE){
		  all(sapply(object, checkOrder, verbose=verbose))
	  })

setMethod("order", signature(...="TrioSetList"),
	  function(..., na.last=TRUE,decreasing=FALSE){
		  x <- list(...)[[1]]
		  for(i in seq_along(x)){
			  x[[i]] <- chromosomePositionOrder(x[[i]])
		  }
		  return(x)
	  })

setMethod("varLabels", signature(object="TrioSetList"),
	  function(object) varLabels(phenoData(object)))

setMethod("pData", signature(object="TrioSetList"),
	  function(object) pData(phenoData(object)))

#' @param name character string of a variable name in the phenoData
#' @aliases $,TrioSetList-method
#' @rdname TrioSetList-class
setMethod("$", signature(x="TrioSetList"),
	  function(x, name){
		  eval(substitute(phenoData(x)$NAME_ARG, list(NAME_ARG=name)))
	  })



setMethod("nrow", signature(x="TrioSetList"), function(x) sum(sapply(lrr(x), nrow)))

setReplaceMethod("featureData", signature(object="TrioSetList", value="list"),
		 function(object, value){
			 object@featureDataList <- value
			 object
		 })

setMethod("gcSubtract", signature(object="TrioSetList"),
	  function(object, ...){
		  .Defunct("methods for GC correction have beem moved to the ArrayTV package available from GitHub")
##		  res <- list()
##		  for(j in seq_along(object)){
##			  res[[j]] <- gcSubtract(object[[j]], ...)
##		  }
##		  object <- stack(RangedDataList(object))
##		  return(object)
	  })

.clone_TrioSetList <- function(object, ids, prefix="clone"){
	if(missing(ids)) ids <- sampleNames(object)
	index <- match(ids, sampleNames(object))
	ids <- as.character(ids)
	r <- lrr(object)
	b <- baf(object)
	rcopy.list <- list()
	bcopy.list <- list()
	for(i in seq_along(r)){
		x <- r[[i]]
		y <- b[[i]]
		open(x)
		open(y)
		rcopy <- initializeBigArray(paste(prefix, "lrr", sep="-"), dim=c(nrow(x), length(ids), 3), vmode="integer")
		bcopy <- initializeBigArray(paste(prefix, "baf", sep="-"), dim=c(nrow(x), length(ids), 3), vmode="integer")
		dimnames(rcopy) <- list(rownames(x),
					colnames(x)[index],
					dimnames(x)[[3]])
		dimnames(bcopy) <- dimnames(rcopy)
		J <- match(ids, colnames(x))
		for(j in seq_along(J)){
			k <- J[j]
			rcopy[, j, ] <- x[, k, ]
			bcopy[, j, ] <- y[, k, ]
		}
		rcopy.list[[i]] <- rcopy
		bcopy.list[[i]] <- bcopy
		close(x)
		close(y)
	}
	adl <- AssayDataList(BAF=bcopy.list, logRRatio=rcopy.list)
	k <- index
	pd <- phenoData(object)[k, ]
	fatherdata <- fatherPhenoData(object)[k, ]
	motherdata <- motherPhenoData(object)[k, ]
	new("TrioSetList",
	    assayDataList=adl,
	    featureDataList=featureData(object),
	    phenoData=pd,
	    fatherPhenoData=fatherdata,
	    motherPhenoData=motherdata,
	    chromosome=chromosome(object),
	    annotation=annotation(object),
	    genome=genomeBuild(object),
	    pedigree=pedigree(object)[k, ])
}


#' @param ranges a \code{GRanges} object
#' @param id a character vector of trio identifiers
#' @param TAUP length-one numeric vector.  Larger values decrease the
#' probability of transitioning to an different state.
#' @param tauMAX the maximum allowed transition probability
#' @param cnStates a length-six numeric vector profiving initial
#' values for the mean copy number for each of the 6 states
#' @param pr.nonmendelian a length-one numeric vector indicating the
#' probability of a non-Mendelian copy number alteration in the offspring
#' @param mdThr a length-one numeric vector indicating the minimum
#' value of the mean minimum distance. Segments with absolute mean
#' value less than \code{mdThr} are not called.
#' @aliases MAP,TrioSetList,GRanges-method
#' @rdname TrioSetList-class
setMethod(MAP, c("TrioSetList", "GRanges"), function(object,
						     ranges,
						     id,
						     TAUP=1e10,
						     tauMAX,
						     cnStates=c(-2, -0.4, 0, 0, 0.4, 1),
						     pr.nonmendelian=1.5e-6,
						     mdThr=0.9, ...){
  ##collapseRanges=TRUE,...){
  .Deprecated("MAP2", msg="This function is deprecated. See MAP2 instead.")
##	.map_trioSetList(object=object,
##			 ranges=ranges,
##			 id=id,
##			 TAUP=TAUP,
##			 tauMAX=tauMAX,
##			 cnStates=cnStates,
##			 pr.nonmendelian=pr.nonmendelian,
##			 mdThr=mdThr, ...)
})

##.map_trioSetList <- function(object,
##			     ranges, id,
##			     TAUP=1e10,
##			     tauMAX,
##			     cnStates=c(-2, -0.4, 0, 0, 0.4, 1),
##			     pr.nonmendelian=1.5e-6,
##			     mdThr=0.9,...){
##  pkgs <- c("VanillaICE", "oligoClasses", "matrixStats", "MinimumDistance")
##  if(missing(id)) id <- sampleNames(object)
##  index.trios <- match(id, sampleNames(object))
##  if(!all(sampleNames(ranges) %in% id))
##    ranges <- ranges[sampleNames(ranges) %in% id, ]
##  if(!all(id %in% sampleNames(ranges))){
##    object <- object[, match(unique(sampleNames(ranges)), id)]
##    id <- id[id %in% sampleNames(ranges)]
##  }
##  chrom.ranges <- unique(chromosome(ranges))
##  chrom.object <- paste("chr", chromosome(object), sep="")
##  object <- object[chrom.object %in% chrom.ranges]
##  ranges <- ranges[chrom.ranges %in% chrom.object, ]
##  ## only call segs that are "nonzero"
##  if("mindist.mad" %in% colnames(elementMetadata(ranges))){
##    mads <- pmax(elementMetadata(ranges)$mindist.mad, .1)
##    abs.thr <- abs(elementMetadata(ranges)$seg.mean)/mads > mdThr
##  } else{
##    ## call all segments
##    abs.thr <- rep(TRUE, length(ranges))
##  }
##  elementMetadata(ranges)$exceeds.md.thr <- abs.thr
##  ocSamples(1) ## has to be 1. This will process 3 samples per alotted CPU
##  chunks <- splitIndicesByLength(index.trios, ocSamples())
##  rlist <- lrr(object)
##  blist <- baf(object)
##  pos <- unlist(position(object))
##  chr <- rep(chromosome(object), elementNROWS(object))
##  build <- genomeBuild(object)
##  sl <- setSequenceLengths(build,
##                           paste("chr", chromosome(object), sep=""))
##  feature.granges <- GRanges(paste("chr", chr, sep=""), IRanges(pos, pos),
##                             seqlengths=sl)
##  grFun <- generatorTransitionProbs(chr, pos, build, TAUP=TAUP, tauMAX=tauMAX)
##  is.snp <- unlist(lapply(featureDataList(object), isSnp))
##  snp.index <- which(is.snp)
##  anyNP <- any(!is.snp)
##  center <- TRUE
##  pkgs <- c("oligoClasses", "VanillaICE")
##  isff <- is(rlist[[1]], "ff")
##  if(isff) pkgs <- c("ff", pkgs)
##  matrixFun <- generatorMatrix(rlist, blist, chr, center=TRUE,
##                               snp.index=snp.index, anyNP=anyNP,
##                               ped=pedigree(object))
##  overlapFun <- generatorOverlapFeatures(feature.granges)
##  grl <- split(ranges, sampleNames(ranges))
##  grl <- grl[match(sampleNames(object), names(grl))]
##  rm(pos, chr, blist, rlist); gc()
##  i <- NULL
##  results <- foreach(i=chunks, granges=grl, .packages=pkgs) %dopar% {
##    emit <- viterbi2Wrapper(index.samples=i,
##                            snp.index=snp.index,
##                            anyNP=anyNP,
##                            is.log=TRUE,
##                            limits=c(-4, 3),
##                            cnStates=cnStates,
##                            grFun=grFun,
##                            matrixFun=matrixFun,
##                            returnEmission=TRUE, ...)
##    granges <- sort(granges)
##    ranges <- loglik(emit=emit,
##                     ranges=granges,
##                     pr.nonmendelian=pr.nonmendelian,
##                     overlapFun=overlapFun)
##    chr.arm <- .getArm(chromosome(ranges), start(ranges), build)
##    ranges <- combineRangesByFactor(ranges, paste(chr.arm, state(ranges), sep="_"))
##    ranges
##  }
##  results <- unlist(GRangesList(results))
##  metadata(results) <- metadata(ranges)
##  return(results)
##}

#' @param md a list of minimum distance matrices. Length of list
#' should be the same as the length of the \code{TrioSetList} object.
#' @param segmentParents logical. Whether to segment the parental log R ratios.
#' @param verbose logical. Whether to display messages indicating progress.
#' @param genome a character vector indicating the UCSC genome build
#' used for the annotation (i.e., 'hg18' or 'hg19').
#' @aliases segment2,TrioSetList-method
#' @rdname TrioSetList-class
setMethod("segment2", signature(object="TrioSetList"),
	  function(object, md=NULL, segmentParents=TRUE, verbose=TRUE, ...){
            segmentTrioSetList(object, md, segmentParents=segmentParents, verbose=verbose, ...)
	  })


#' @param pos a list of the genomic positions (integers)
#' @param chrom list of chromosome names
#' @param featureNames a list of the marker names
#' @aliases segment2,list-method
#' @rdname TrioSetList-class
setMethod("segment2", signature(object="list"),
	  function(object, pos, chrom, id=NULL, featureNames, segmentParents=TRUE, verbose=TRUE, genome, ...){
            ## elements of list must be a matrix or an array
            if(missing(genome)) stop("must specify UCSC genome build")
            segs <- segmentList(object, pos, chrom, id, featureNames, segmentParents=segmentParents, verbose=verbose, genome=genome, ...)
            metadata(segs) <- list(genome=genome)
            segs
	  })
