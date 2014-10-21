setMethod("initialize", "TrioSet",
	  function(.Object,
		   assayData=assayDataNew(logRRatio=logRRatio, BAF=BAF, ...),
		   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
		   fatherPhenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
		   motherPhenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
		   annotation=character(),
		   featureData=GenomeAnnotatedDataFrameFrom(assayData, annotation, genome=genome),
		   experimentData=new("MIAME"),
		   protocolData=phenoData[, integer(0)],
		   logRRatio=array(NA, dim=c(0, 0, 3)),
		   BAF=array(NA, dim=c(0,0,3)),
		   pedigree=Pedigree(),
		   mindist=NULL,
		   genome=c("hg19", "hg18"), ...){
		  .Object@pedigree <- pedigree
		  .Object@fatherPhenoData <- fatherPhenoData
		  .Object@motherPhenoData <- motherPhenoData
		  callNextMethod(.Object,
				 assayData=assayData,
				 phenoData=phenoData,
				 ##fatherPhenoData=fatherPhenoData,
				 motherPhenoData=motherPhenoData,
				 featureData=featureData,
				 experimentData=experimentData,
				 annotation=annotation,
				 protocolData=protocolData,
				 pedigree=pedigree,
				 mindist=mindist,
				 genome=match.arg(genome), ...)
	  })

setValidity("TrioSet", function(object){
	ped <- pedigree(object)
	validObject(ped)
	validObject(featureData(object))
	nms <- ls(assayData(object))
	if(!all(c("BAF", "logRRatio") %in% nms)){
		msg <- "BAF and logRRatio are required elements of the assayData"
		return(msg)
	}
	elt <- nms[[1]]
	elt <- assayData(object)[[elt]]
	if(ncol(elt) > 0){
		sns.ped <- sampleNames(ped)
		if(length(sns.ped) != ncol(elt)){
			return("Number of samples in pedigree slot should be the same as the number of columns in the TrioSet object")
		}
	}
	if(!identical(sampleNames(object), sampleNames(phenoData(object)))){
		return("sampleNames of TrioSetList object must be the same as the sampleNames of the phenoData")
	}
	if(!identical(fatherNames(object), originalNames(sampleNames(fatherPhenoData(object))))){
		return("fatherNames of TrioSetList object must be the same as the sampleNames of the fatherPhenoData")
	}
	if(!identical(motherNames(object), originalNames(sampleNames(motherPhenoData(object))))){
		stop("motherNames of TrioSetList object must be the same as the sampleNames of the motherPhenoData")
	}
	if(!is.null(mindist(object))){
		if(!identical(colnames(mindist(object)), sampleNames(object)))
			stop("colnames of mindist matrix must be same as the sampleNames of the TrioSet object")
	}
})

setMethod("updateObject", signature(object="TrioSet"),
	  function(object, ..., verbose=FALSE){
		  if (verbose) message("updateObject(object = 'TrioSetList')")
		  if(!is(featureData(object), "GenomeAnnotatedDataFrame")){
			  featureData(object) <- updateObject(featureData(object))
		  }
		  return(object)
	  })


#' @param object a \code{TrioSet} object
#' @aliases pedigree,TrioSet-method
#' @rdname TrioSet-class
setMethod("pedigree", signature(object="TrioSet"), function(object) object@pedigree)


setMethod("lrr", "TrioSet", function(object) assayDataElement(object, "logRRatio"))
setReplaceMethod("lrr", c("TrioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "logRRatio", value)
	 })

setMethod("baf", "TrioSet",
	  function(object) {
		  assayDataElement(object, "BAF")
	 })
setReplaceMethod("baf", c("TrioSet", "array"),
		 function(object, value) {
			 assayDataElementReplace(object, "BAF", value)
	 })
setReplaceMethod("baf", c("TrioSet", "ff_array"),
		 function(object, value) {
			 assayDataElementReplace(object, "BAF", value)
	 })

setMethod("fatherPhenoData", signature(object="TrioSet"),
	  function(object) object@fatherPhenoData)
setMethod("motherPhenoData", signature(object="TrioSet"),
	  function(object) object@motherPhenoData)
setMethod("offspringPhenoData", signature(object="TrioSet"),
	  function(object) phenoData(object))

#' Deprecated constructor for \code{TrioSet} class
#'
#' The \code{TrioSet} class has been deprecated and may be removed in
#' a future release.
#'
#' @param pedigreeData an object of class \code{Pedigree}
#' @param sample.sheet a \code{data.frame} containing metadata on the trios
#' @param row.names a character vector providing row identifiers for
#' the \code{sample.sheet} argument that match the names of the
#' trios in the \code{pedigreeData} argument.
#' @param lrr a matrix of log R ratios
#' @param baf a matrix of B allele frequencies
#' @param featureData a \code{GenomeAnnotatedDataFrame} object for the SNPs/nonpolymorphic markers
#' @param cdfname character string indicating the annotation package used to extract physical position and chromosome of markers
#' @param drop logical.  When FALSE, the dimnames on the log R ratio and BAF arrays is set to NULL
#' @param mindist  can be either NULL or a matrix of the minimum distance
#' @param genome character string providing the UCSC genome build
#' @return \code{TrioSet}
#' @examples
#' 	path <- system.file("extdata", package="MinimumDistance")
#' 	load(file.path(path, "logRratio.rda"))
#' 	load(file.path(path, "baf.rda"))
#' 	load(file.path(path, "pedigreeInfo.rda"))
#' 	trioSet <- TrioSet(lrr=logRratio,
#' 			   baf=baf,
#' 			   pedigree=Pedigree(pedigreeInfo),
#' 			   cdfname="human610quadv1bCrlmm",
#' 			   genome="hg18")
#' @export
TrioSet <- function(pedigreeData=Pedigree(),
		    sample.sheet,
		    row.names=NULL,
		    lrr,
		    baf,
		    featureData,
		    cdfname,
		    drop=TRUE,
		    mindist=NULL, genome=c("hg19", "hg18")){
  if(missing(lrr) | missing(baf)){
    object <- new("TrioSet",
                  pedigree=pedigreeData)
    return(object)
  } else{
    if(ncol(lrr) > 0 & nrow(pedigreeData)==0)
      stop("pedigreeData has zero rows")
  }
  if(!missing(lrr) & !missing(baf)){
    if(!identical(rownames(lrr), rownames(baf)))
      stop("rownames of lrr and baf are not identical")
    if(!identical(dim(lrr), dim(baf)))
      stop("lrr and baf must have the same dimension")
    if(!(is(lrr[1,1], "integer") & is(baf[1,1], "integer"))){
      stop("rr and baf must be integers. Use integerMatrix(x, scale=100) to transform log R ratios and integerMatrix(x, scale=1000) for B allele frequencies")
    }
  }
  if(missing(featureData)){
    if(missing(cdfname)) stop("If featureData is not supplied, a valid cdfname must be provided for feature annotation")
    featureData <- GenomeAnnotatedDataFrameFrom(lrr, cdfname, genome=match.arg(genome))
    fD <- featureData[order(chromosome(featureData), position(featureData)), ]
    rm(featureData); gc()
  } else {
    if(!is(featureData, "AnnotatedDataFrame")) stop("featureData must be an AnnotatedDataFrame or a GenomeAnnotatedDataFrame")
    fD <- featureData
  }
  is.present <- sampleNames(fD) %in% rownames(lrr)
  if(!all(is.present)) fD <- fD[is.present, ]
  if(!is.null(rownames(lrr))){
    index <- match(sampleNames(fD), rownames(lrr))
    if(length(index) == 0) {
      if(!missing(cdfname)){
        msg <- paste("rownames for log R ratios do not match feature ids with annotation package ", cdfname)
        stop(msg)
      }
    }
    lrr <- lrr[index, ]
    baf <- baf[index, ]
    stopifnot(all(identical(rownames(lrr), sampleNames(fD))))
  }
  np <- nrow(trios(pedigreeData))
  trio.names <- array(NA, dim=c(length(offspringNames(pedigreeData)), 1, 3))
  dimnames(trio.names) <- list(offspringNames(pedigreeData), "sampleNames", c("F", "M", "O"))
  trio.names[, "sampleNames", ] <- as.matrix(trios(pedigreeData))
  father.names <- fatherNames(pedigreeData)
  mother.names <- motherNames(pedigreeData)
  offspring.names <- offspringNames(pedigreeData)
  father.index <- match(father.names, colnames(lrr))
  if(length(father.index) == 0) stop("father ids in pedigree do not match any of the column names of the lrr matrix")
  mother.index <- match(mother.names, colnames(lrr))
  if(length(mother.index) == 0) stop("mother ids in pedigree do not match any of the column names of the lrr matrix")
  offspring.index <- match(offspring.names, colnames(lrr))
  if(length(offspring.index) == 0) stop("offspring ids in pedigree do not match any of the column names of the lrr matrix")
  nr <- nrow(lrr)
  np <- length(offspring.names)
  bafArray <- initializeBigArray("baf", dim=c(nr, np, 3), vmode="integer")
  logRArray <- initializeBigArray("lrr", dim=c(nr, np, 3), vmode="integer")
  dimnames(bafArray)[[3]] <- dimnames(logRArray)[[3]] <- c("F", "M", "O")
  logRArray[,,"F"] <- lrr[, father.index]
  logRArray[,,"M"] <- lrr[, mother.index]
  logRArray[,,"O"] <- lrr[, offspring.index]
  bafArray[,,"F"] <- baf[, father.index]
  bafArray[,,"M"] <- baf[, mother.index]
  bafArray[,,"O"] <- baf[, offspring.index]
  if(!drop){
    dimnames(bafArray)[c(1,2)] <- dimnames(logRArray)[c(1,2)] <- list(sampleNames(fD), colnames(lrr)[offspring.index])
  }
  if(nrow(pedigreeData) > 0){
    if(!missing(sample.sheet)){
      if(is.null(row.names)){
        row.names <- rownames(sample.sheet)
      }
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
  object <- new("TrioSet",
                BAF=bafArray,
                logRRatio=logRArray,
                phenoData=phenoData,
                fatherPhenoData=fatherPhenoData,
                motherPhenoData=motherPhenoData,
                pedigree=pedigreeData,
                featureData=fD,
                mindist=mindist,
                genome=match.arg(genome))
}


#' @aliases show,TrioSet-method
#' @rdname TrioSet-class
setMethod("show", signature(object="TrioSet"),
	  function(object){
              cat(class( object ), " (storageMode: ", storageMode(object), ")\n", sep="")
              adim <- dim(object)
	      cat("assayData:\n")
              cat("  element names:",
		  paste(assayDataElementNames(object), collapse=", "), "\n")
	      cat("  dimension:\n")
	      print(adim)
	      cat("genome:", genomeBuild(object), "\n")
	  })

setMethod("open", "TrioSet", function(con, ...){
	object <- con
	if(!isFF(object)) return()
	names <- ls(assayData(object))
	L <- length(names)
	for(i in 1:L) open(eval(substitute(assayData(object)[[NAME]], list(NAME=names[i]))))
	L <- length(names)
	if("MAD" %in% varLabels(object)){
		if(is(object$MAD, "ff")) open(object$MAD)
	}
	open(mindist(object))
	return(TRUE)
})

setMethod("close", "TrioSet", function(con, ...){
	##browser()
	##con is just to keep the same generic arguments
	object <- con
	if(!isFF(object)) return()
	names <- ls(assayData(object))
	L <- length(names)
	for(i in 1:L) close(eval(substitute(assayData(object)[[NAME]], list(NAME=names[i]))))
	if("MAD" %in% varLabels(object)){
		if(is(object$MAD, "ff")) close(object$MAD)
	}
	close(mindist(object))
	return()
})


setReplaceMethod("sampleNames", signature(object="TrioSet"), function(object, value){
	callNextMethod(object, value)
})

#' @aliases mindist,TrioSet-method
#' @rdname TrioSet-class
setMethod("mindist", "TrioSet", function(object) object@mindist)

#' @param value a \code{matrix}
#' @aliases mindist<-,TrioSet,matrix-method
#' @rdname TrioSet-class
setReplaceMethod("mindist", signature(object="TrioSet", value="matrix"),
		 function(object, value){
			 object@mindist <- value
			 return(object)
		 })

#' @param x a \code{TrioSet} object
#' @aliases dim,TrioSet-method
#' @rdname TrioSet-class
setMethod("dim", "TrioSet", function(x) {
  adim <- callNextMethod(x)
  names(adim) <- c("Features", "Trios", "Members")
  adim
})

setMethod("ncol", signature(x="TrioSet"), function(x) dim(x)[[2]])

#' @aliases trios,TrioSet-method
#' @rdname TrioSet-class
setMethod("trios", signature(object="TrioSet"),
	  function(object){
		  trios(pedigree(object))
                })

#' @param i a numeric vector for subsetting rows  (optional)
#' @param j a numeric vector for subsetting trios (optional)
#' @param ... additional arguments passed to subsetting methods for matrices and data frames
#' @param drop logical. Whether to simplify matrices to numeric
#' vectors.  This should be left as FALSE.
#' @aliases "[",TrioSet,ANY-method
#' @rdname TrioSet-class
setMethod("[", "TrioSet", function(x, i, j, ..., drop = FALSE) {
	if (missing(drop))
		drop <- FALSE
	##	if(length(list(...)) > 0){
	##		k <- list(...)[[1]]
	##	} else k <- NULL
	##	if (missing(i) && missing(j) && is.null(k)) {
	if(missing(i) && missing(j)){
		return(x)
		##if (length(list(...))!=0)
		##	stop("specify features, trios, or samples to subset; use '",
		##	     substitute(x), "$", names(list(...))[[1]],
		##	     "' to access phenoData variables")
		##return(x)
	}
	if (!missing(j) & missing(i)) {
		phenoData(x) <- phenoData(x)[j,, ..., drop = drop]
		protocolData(x) <- protocolData(x)[j,, ..., drop = drop]
		##x@sampleSheet <- x@sampleSheet[j, , , drop=drop]
		##tmp <- pedigree(x)[j, , drop=drop]
		x@pedigree <- pedigree(x)[j, , drop=drop]
		b <- baf(x)[, j, , drop=drop]
		r <- lrr(x)[, j, , drop=drop]
		x <- assayDataElementReplace(x, "logRRatio", r)
		x <- assayDataElementReplace(x, "BAF", b)
		x@fatherPhenoData <- fatherPhenoData(x)[j, ]
		x@motherPhenoData <- motherPhenoData(x)[j, ]
		if(!is.null(mindist(x))){
			mindist(x) <- mindist(x)[, j, drop=FALSE]
		}
		##mad.sample(x) <- mad.sample(x)[j,,...,drop=drop]
	}
	if (!missing(i) & !missing(j)){
		phenoData(x) <- phenoData(x)[j, ..., drop=drop]
		protocolData(x) <- protocolData(x)[j, ..., drop=drop]
		featureData(x) <- featureData(x)[i, ,..., drop=drop]
		b <- baf(x)[i, j, , drop=drop]
		r <- lrr(x)[i, j, , drop=drop]
		##x@sampleSheet <- x@sampleSheet[j, , , drop=drop]
		x@pedigree <- x@pedigree[j, , drop=drop]
		x <- assayDataElementReplace(x, "logRRatio", r)
		x <- assayDataElementReplace(x, "BAF", b)
		x@fatherPhenoData <- fatherPhenoData(x)[j, ]
		x@motherPhenoData <- motherPhenoData(x)[j, ]
		if(!is.null(mindist(x))){
			mindist(x) <- mindist(x)[i, j, drop=FALSE]
		}
	}
	if(!missing(i) & missing(j)){
		featureData(x) <- featureData(x)[i, ,..., drop=drop]
		b <- baf(x)[i, , , drop=drop]
		r <- lrr(x)[i, , , drop=drop]
		x <- assayDataElementReplace(x, "logRRatio", r)
		x <- assayDataElementReplace(x, "BAF", b)
		if(!is.null(mindist(x))){
			mindist(x) <- mindist(x)[i, , drop=FALSE]
		}
	}
	return(x)
})

setMethod("checkOrder", signature(object="TrioSet"),
	  function(object, verbose=FALSE){
            .checkOrder(object, verbose)
	  })

setMethod("todf", signature(object="TrioSet", rangeData="RangedData"),
	  function(object, rangeData, frame, ...){
		  ## FIX
		  stop("requires mindist(object)")
		  stopifnot(nrow(rangeData) == 1)
		  if(missing(frame)){
			  w <- width(rangeData)
			  frame <- w/0.05  * 1/2
		  }
		  ##overlaps <- findOverlaps(object, rangeData, max.gap=frame)
		  overlaps <- findOverlaps(rangeData, featureData(object), maxgap=frame)
		  marker.index <- subjectHits(overlaps)
		  ##marker.index <- featuresInRangeData(object, rangeData, FRAME=frame)
		  id <- rangeData$id
		  sample.index <- match(id, sampleNames(object))
		  stopifnot(length(sample.index)==1)
		  is.ff <- is(lrr(object), "ff")
		  if(is.ff){
			  open(baf(object))
			  open(lrr(object))
			  open(mindist(object))
		  }
		  b <- baf(object)[marker.index, sample.index, ]
		  r <- lrr(object)[marker.index, sample.index, ]
		  md <- mindist(object)[marker.index, sample.index]
		  if(is.ff){
			  close(baf(object))
			  close(lrr(object))
			  close(mindist(object))
		  }
		  id <- matrix(c("father", "mother", "offspring"), nrow(b), ncol(b), byrow=TRUE)
		  empty <- rep(NA, length(md))
		  ## A trick to add an extra panel for genes and cnv
		  ##df <- rbind(df, list(as.integer(NA), as.numeric(NA), as.numeric(NA), as.factor("genes")))
		  ## The NA's are to create extra panels (when needed for lattice plotting)
		  id <- c(as.character(id), rep("min dist",length(md)))##, c("genes", "CNV"))
		  b <- c(as.numeric(b), empty)
		  r <- c(as.numeric(r), md)
		  x <- rep(position(object)[marker.index], 4)/1e6
		  is.snp <- rep(isSnp(object)[marker.index], 4)
		  df <- data.frame(x=x, b=b, r=r, id=id, is.snp=is.snp)
		  df2 <- data.frame(id=c(as.character(df$id), "genes", "CNV"),
				    b=c(df$b, NA, NA),
				    r=c(df$r, NA, NA),
				    x=c(df$x, NA, NA),
				    is.snp=c(df$is.snp,NA, NA))
		  df2$id <- factor(df2$id, levels=c("father", "mother", "offspring", "min dist", "genes", "CNV"), ordered=TRUE)
		  return(df2)
	  })

setMethod("prune", signature(object="TrioSet", ranges="RangedDataCNV"),
	  function(object, ranges, md, verbose=TRUE, ...){
		  pruneTrioSet(object=object, ranges=ranges, md=md, verbose=verbose, ...)
	 })

pruneTrioSet <- function(object, ranges, md, verbose=TRUE, ...){
	CHR <- unique(chromosome(object))
	if(verbose) message("Pruning chromosome ", CHR)
	id <- unique(sampleNames(ranges))
	index <- which(chromosome(ranges) == CHR & sampleNames(ranges) %in% id)
	ranges <- ranges[index, ]
	rdList <- vector("list", length(unique(id)))
	if(verbose){
		message("\tPruning ", length(unique(id)), " files.")
		pb <- txtProgressBar(min=0, max=length(unique(id)), style=3)
	}
	for(j in seq_along(id)){
		if (verbose) setTxtProgressBar(pb, j)
		sampleId <- id[j]
		rd <- ranges[sampleNames(ranges) == sampleId, ]
		stopifnot(nrow(rd) > 0)
		## assign the mad of the minimum distance to the ranges
		k <- match(sampleId, sampleNames(object))
		##rd$mad <- object[[1]]$mindist.mad[k]
		genomdat <- md[, k]
		##genomdat <- as.numeric(mindist(object)[, k])/100
		## This function currently returns a RangedData object
		rdList[[j]] <- pruneMD(genomdat, rd,  ...)
	}
	if(verbose) close(pb)
	rd <- stack(RangedDataList(rdList))
	rd <- rd[, -ncol(rd)]
	return(rd)
}

##setMethod("offspringNames", signature(object="TrioSet"), function(object){
##	phenoData2(object)[, "id", "O"]
##})
##setReplaceMethod("offspringNames", signature(object="TrioSet", value="character"),
##		 function(object, value){
##			 phenoData2(object)[, "id", "O"] <- value
##			 object
##		 })

setMethod("fatherNames", signature(object="TrioSet"), function(object){
	##phenoData2(object)[, "id", "F"]
	fatherNames(pedigree(object))
})
##setReplaceMethod("fatherNames", signature(object="TrioSet", value="character"),
##		 function(object, value){
##			 phenoData2(object)[, "id", "F"] <- value
##			 object
##		 })
setMethod("motherNames", signature(object="TrioSet"), function(object){
	##phenoData2(object)[, "id", "M"]
	motherNames(pedigree(object))
})
##setReplaceMethod("motherNames", signature(object="TrioSet", value="character"),
##		 function(object, value){
##			 phenoData2(object)[, "id", "M"] <- value
##			 object
##		 })
##fmoNames <- function(object){
##	tmp <- cbind(fatherNames(object), motherNames(object), offspringNames(object))
##	colnames(tmp) <- c("F", "M", "O")
##	return(tmp)
##}



# setMethod("xyplot", signature(x="formula", data="TrioSet"),
# 	  function(x, data, ...){
#             if("range" %in% names(list(...))){
#               res <- xyplot2(x, data, ...)
#             } else {
#               callNextMethod()
#             }
# 	  })


##setMethod("phenoData2", signature(object="TrioSet"),
##	  function(object) object@phenoData2)
setMethod("allNames", signature(object="TrioSet"), function(object) allNames(pedigree(object)))


setMethod("order", "TrioSet", ##signature(...="TrioSet"),
	  function(..., na.last=TRUE, decreasing=FALSE){
            x <- list(...)[[1]]
            chromosomePositionOrder(x)
	  })


#' @param verbose logical. Whether to display messages indicating progress.
#' @aliases calculateMindist,TrioSet-method
#' @rdname calculateMindist
setMethod("calculateMindist", signature(object="TrioSet"),
	  function(object, verbose=TRUE, ...){
		  calculateMindist(lrr(object))
	  })

setMethod("gcSubtract", signature(object="TrioSet"),
	  function(object, method=c("speed", "lowess"), trio.index, ...){
		  .Defunct("methods for GC correction have been moved to the ArrayTV package available from GitHub")
	  })


#' @param ranges a \code{GRanges} object
#' @param transition_param an object of class \code{TransitionParam}
#' @param emission_param an object of class \code{EmissionParam}
#' @param mdThr the minimum absolute value of the minimum distance
#' segment mean. Segments with means below \code{mdThr} in absolute
#' value will not be called as they are unlikely to be de novo.
#' @aliases MAP,TrioSet,GRanges-method
#' @rdname TrioSet-class
setMethod(MAP, c("TrioSet", "GRanges"), function(object,
						 ranges,
                                                 ##id,
                                                 transition_param=TransitionParam(),
                                                 emission_param=EmissionParam(),
						 mdThr=0.9, ...){
  .Deprecated("MAP2", msg="This function is deprecated. See MAP2 instead.")
})





#' @param md a matrix of the minimum distance
#' @param segmentParents logical.  Whether to segment the log R ratios
#' of the parents using circular binary segmentation.
#' @param verbose logical. Whether to display messages that indicate progress.
#' @aliases segment2,TrioSet-method
#' @seealso \code{\link[DNAcopy]{segment}}
#' @rdname segment2
setMethod("segment2", signature(object="TrioSet"),
	  function(object, md=NULL, segmentParents=TRUE, verbose=TRUE, ...){
            segmentTrioSet(object, md=md, segmentParents=segmentParents, verbose=verbose, ...)
	  })


#' @aliases segment2,matrix-method
#' @rdname segment2
setMethod("segment2", signature(object="matrix"),
	  function(object, pos, chrom, id, featureNames, ...){
            stopifnot(is(id, "character"))
            segmentMatrix(object, pos, chrom, id, featureNames, ...)
	  })

#' @aliases segment2,ff_matrix-method
#' @rdname segment2
setMethod("segment2", signature(object="ff_matrix"),
	  function(object, pos, chrom, id, featureNames, ...){
            segmentff_matrix(object, pos, chrom, id, featureNames, ...)
            ##segs <- foreach(i=seq_along(ilist), .packages="MinimumDistance") %dopar% segmentMatrix(object[, ilist[[i]]], pos=pos, chrom=chrom, id=id[ilist[[i]]], featureNames, ...)
	  })

#' @param featureNames character vector specifying marker names for subsetting \code{object}
#' @param id character vector of trio identifiers for subsetting \code{object}
#' @param chrom character or integer vector of chromosome names
#' @param pos integer vector of physical position of markers in the genome
#' @aliases segment2,arrayORff_array-method
#' @rdname segment2
setMethod("segment2", signature(object="arrayORff_array"),
	  function(object, pos, chrom, id, featureNames, segmentParents=TRUE, verbose=TRUE, ...){
            segmentArray(object, pos, chrom, id, featureNames, segmentParents=segmentParents, verbose=verbose, ...)
	  })
