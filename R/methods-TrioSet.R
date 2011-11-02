setMethod("initialize", signature(.Object="TrioSet"),
	  function(.Object,
		   logRRatio=array(NA, dim=c(0L,0L,1L)), ## Need the one so that annotatedDataFrameFromArray will work...
		   BAF=array(NA, dim=dim(logRRatio)),
		   phenoArray=array(NA, dim=dim(logRRatio)), ##perhaps replace with initializeBigArray
		   mindist=initializeBigMatrix("mindist", nr=nrow(logRRatio), ncol(logRRatio), vmode="double"),##perhaps replace with initializeBigMatrix
		   mad=matrix(NA, nrow(logRRatio), ncol(logRRatio)),
		   ...){
		  .Object@phenoData2 <- phenoArray
		  .Object <- callNextMethod(.Object, logRRatio=logRRatio, BAF=BAF, ...)
		  ##.Object <- callNextMethod(.Object, logRRatio=logRRatio, BAF=BAF, ...)
		  ##.Object <- callNextMethod(.Object, BAF=BAF, ...)
		  ##.Object <- callNextMethod(.Object)
		  ##.Object <- assayDataElementReplace(.Object, "logRRatio", array(NA, dim=c(0,0,0)))
		  .Object@mindist <- mindist
		  .Object@mad <- mad
		  return(.Object)
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
	object <- callNextMethod(object, value)
	if(!is.null(mindist(object))){
		colnames(mindist(object)) <- value
	}
	if(!is.null(mad(object))){
		##colnames(mad(object)) <- value
		rownames(mad(object)) <- value
	}
	return(object)
})

setMethod("show", signature(object="TrioSet"), function(object){
	cat(class( object ), " (storageMode: ", storageMode(object), ")\n", sep="")
	adim <- dim(object)
	if (length(adim)>1)
		cat("assayData:",
		    if (length(adim)>1)
		    paste(adim[[1]], "Features,",
			  adim[[2]], "Trios, ",
			  adim[[3]], "Family members (Father, Mother, Offspring)") else NULL,
		    "\n")
	cat("  element names:",
	    paste(assayDataElementNames(object), collapse=", "), "\n")
	Biobase:::.showAnnotatedDataFrame(protocolData(object),
				labels=list(object="protocolData"))
	Biobase:::.showAnnotatedDataFrame(phenoData(object),
				labels=list(object="phenoData"))
	Biobase:::.showAnnotatedDataFrame(featureData(object),
				labels=list(
				object="featureData",
				sampleNames="featureNames",
				varLabels="fvarLabels",
				varMetadata="fvarMetadata"))
	cat("experimentData: use 'experimentData(object)'\n")
	pmids <- pubMedIds(object)
	if (length(pmids) > 0 && all(pmids != ""))
		cat("  pubMedIds:", paste(pmids, sep=", "), "\n")
	cat("Annotation:", annotation(object), "\n")
	cat("mindist:", class(mindist(object)), "\n")
	cat("phenoData2:", class(object@phenoData2), "\n")
	cat("mad:", class(mad(object)), "\n")
})

setReplaceMethod("sampleNames", signature(object="TrioSet"),
		 function(object,value){
			 object <- callNextMethod(object, value)
			 row.names(object@phenoData2) <- value
			 object
		 })

setMethod("mindist", "TrioSet", function(object) object@mindist)
setReplaceMethod("mindist", signature(object="TrioSet", value="ANY"),
		 function(object, value){
			 object@mindist <- value
			 return(object)
		 })

##setReplaceMethod("trioNames", signature(object="TrioSet"),
##		 function(object,value){
##			 object <- callNextMethod(object, value)
##			 row.names(object@phenoData2) <- value
##			 object
##		 })

setMethod("dim", "TrioSet", function(x) assayDataDim(assayData(x)))


## Fix strange behavor for [ with ff_arrays.
setMethod("[", signature(x="ff_array"), function(x, i, j, ..., drop=FALSE){
	if(missing(drop)) drop <- TRUE
	if(length(list(...)) > 0){
		k <- list(...)[[1]]
		if(is(k, "character")) {
			k <- match(dimnames(x)[[3]])
			stopifnot(length(k)>0)
		}
	} else k <- NULL
	if(is.null(k)){
		if(!missing(i) & missing(j)){
			x <- x[i, 1:ncol(x), 1:3, ..., drop=drop]
		}
		if(!missing(i) & !missing(j)){
			x <- x[1:nrow(x), j, 1:3, drop=drop]
		}
		if(missing(i) & !missing(j)){
			x <- x[1:nrow(x), j, 1:3, drop=drop]
		}
	} else {
		if(!missing(i) & missing(j)){
			x <- x[i, 1:ncol(x), k, drop=drop]
		}
		if(!missing(i) & !missing(j)){
			x <- x[i, j, k, drop=drop]
		}
		if(missing(i) & !missing(j)){
			x <- x[1:nrow(x), j, k, drop=drop]
		}
	}
	return(x)
})

setMethod("[", "TrioSet", function(x, i, j, ..., drop = FALSE) {
	if (missing(drop))
		  drop <- FALSE
	if(length(list(...)) > 0){
		k <- list(...)[[1]]
	} else k <- NULL
	if (missing(i) && missing(j) && is.null(k)) {
		if (length(list(...))!=0)
			stop("specify features, trios, or samples to subset; use '",
			     substitute(x), "$", names(list(...))[[1]],
			     "' to access phenoData variables")
		return(x)
	}
	if (!missing(j)) {
		phenoData(x) <- phenoData(x)[j,, ..., drop = drop]
		protocolData(x) <- protocolData(x)[j,, ..., drop = drop]
		mad.sample(x) <- mad.sample(x)[j,,...,drop=drop]
	}
	if (!missing(i))
		featureData(x) <- featureData(x)[i, ,..., drop=drop]
	if (!is.null(k) && !missing(j)){
		x@phenoData2 <- x@phenoData2[j, k, , drop=drop]
	}
	if(!missing(j) && is.null(k)){
		x@phenoData2 <- x@phenoData2[j, , , drop=drop]
	}
	if(missing(j) && !is.null(k)){
		x@phenoData2 <- x@phenoData2[, k, , drop=drop]
	}
	if(!missing(i) & missing(j)){
		mindist(x) <- mindist(x)[i, ,drop=drop]
	}
	if(!missing(i) & !missing(j)){
		mindist(x) <- mindist(x)[i, j, drop=drop]
	}
	if(missing(i) & !missing(j)){
		mindist(x) <- mindist(x)[, j, drop=drop]
	}
	## assayData; implemented here to avoid function call
	orig <- assayData(x)
	storage.mode <- Biobase:::assayDataStorageMode(orig)
	if(is.null(k)){
		assayData(x) <-
			switch(storage.mode,
			       environment =,
			       lockedEnvironment = {
				       aData <- new.env(parent=emptyenv())
				       if (missing(i))                     # j must be present
					       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][, j, , ..., drop = drop]
				       else {                              # j may or may not be present
					       if (missing(j))
						       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i,, , ..., drop = drop]
					       else
						       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i, j, , ..., drop = drop]
				       }
				       if ("lockedEnvironment" == storage.mode) Biobase:::assayDataEnvLock(aData)
				       aData
			       },
			       list = {
				       if (missing(i))                     # j must be present
					       lapply(orig, function(obj) obj[, j, , ..., drop = drop])
				       else {                              # j may or may not be present
					       if (missing(j))
						       lapply(orig, function(obj) obj[i,, , ..., drop = drop])
					       else
						       lapply(orig, function(obj) obj[i, j, , ..., drop = drop])
				       }
			       })
		return(x)
	} else{ ## not missing k
		assayData(x) <-
			switch(storage.mode,
			       environment =,
			       lockedEnvironment = {
				       aData <- new.env(parent=emptyenv())
				       if (missing(i))                     # j must be present
					       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][, j, k, drop = drop]
				       else {                              # j may or may not be present
					       if (missing(j))
						       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i,, k, drop = drop]
					       else
						       for(nm in ls(orig)) aData[[nm]] <- orig[[nm]][i, j, k, drop = drop]
				       }
				       if ("lockedEnvironment" == storage.mode) Biobase:::assayDataEnvLock(aData)
				       aData
			       },
			       list = {
				       if (missing(i))                     # j must be present
					       lapply(orig, function(obj) obj[, j, k, ..., drop = drop])
				       else {                              # j may or may not be present
					       if (missing(j))
						       lapply(orig, function(obj) obj[i,, k, ..., drop = drop])
					       else
						       lapply(orig, function(obj) obj[i, j, k, ..., drop = drop])
				       }
			       })
		return(x)
	}
})



setReplaceMethod("logR", signature(object="TrioSet", value="ANY"),
		 function(object, value){
			 assayDataElementReplace(object, "logRRatio", value)
		 })
setReplaceMethod("baf", signature(object="TrioSet", value="ANY"),
		 function(object, value){
			 assayDataElementReplace(object, "BAF", value)
		 })

fullId <- function(object) object@phenoData2[, "id", ]

setMethod("calculateMindist", signature(object="TrioSet"),
	  function(object, ..., verbose=TRUE){
        sns <- sampleNames(object)
	is.ff <- is(logR(object), "ff")
	if(is.ff){
		invisible(open(lrr(object)))
		if(!is.null(mindist(object)))
			invisible(open(mindist(object)))
	}
	## only instantiate a new object if slot is NULL
	if(is.null(mindist(object))){
		md <- initializeBigMatrix("mindist", nr=nrow(object), nc=ncol(object),
					  vmode="double")
	} else md <- mindist(object)
	if(verbose){
		message("\t\tComputing Bayes factors for ", ncol(object), " files.")
		pb <- txtProgressBar(min=0, max=ncol(object), style=3)
	}
	for(j in seq(length=ncol(object))){
		##if(j %% 100 == 0) cat(".")
		if(verbose) setTxtProgressBar(pb, j)
		lr <- logR(object)[, j, ]
		d1 <- lr[, "F"] - lr[, "O"]
		d2 <- lr[, "M"] - lr[, "O"]
		I <- as.numeric(abs(d1) <= abs(d2))
		md[, j] <- I*d1 + (1-I)*d2
		##mindist(object)[, j] <- md
		##object$MAD[j] <- mad(md, na.rm=TRUE)
	}
	if(verbose) close(pb)
	if(is.ff){
		if(!is.null(mindist(object)))
			close(mindist(object))
		close(lrr(object))
	}
	##return(object)
	dimnames(md) <- list(featureNames(object), sampleNames(object))
	return(md)
})

setMethod("mad", signature(x="TrioSet"), function(x) x@mad)
setMethod("mad.marker", signature(x="TrioSet"), function(x) fData(x)$marker.mad)
setMethod("mad.sample", signature(x="TrioSet"), function(x) x@mad)
setReplaceMethod("mad.sample", signature(x="TrioSet", value="ANY"),
	  function(x, value){
		  x@mad <- value
		  return(x)
	  })
setReplaceMethod("mad.marker", signature(x="TrioSet", value="ANY"),
	  function(x, value){
		  fData(x)$marker.mad <- value
		  return(x)
	  })

setMethod("phenoData2", signature(object="TrioSet"), function(object) object@phenoData2)
setReplaceMethod("phenoData2", signature(object="TrioSet", value="ANY"), function(object, value){
	object@phenoData2 <- value
	object
})
setMethod("varLabels2", signature(object="TrioSet"), function(object) colnames(phenoData2(object)))

offspringForId <- function(id, pedigree){
	## we have a vector of ids.
	index <- match(id, pedigree.char)
}


setMethod("segment2", signature(object="matrix"),
	  function(object, pos, chrom, ...){
		  res <- segmentMatrix(object, pos, chrom, ...)
	  })

setMethod("segment2", signature(object="array", id="data.frame"),
	  function(object, pos, chrom, id, ...){
		  stopifnot(length(dim(object))==3) ## expects a 3d array
		  J <- dim(object)[[3]]
		  resList <- vector("list", J)
		  if("verbose" %in% names(list(...))){
			  verbose <- list(...)[["verbose"]]
		  }
		  for(j in seq_len(J)){
			  if(verbose) message("Processing ", ncol(object), " samples")
			  obj <- object[, , j]
			  colnames(obj) <- id[, j]
			  resList[[j]] <- segmentMatrix(obj, pos, chrom, ...)
			  rm(obj)
		  }
		  res <- stack(RangedDataList(resList))
		  j <- match("sample", colnames(res))
		  if(length(j) == 1) res <- res[, -j]
		  return(res)
	  })

setMethod("segment2", signature(object="list", pos="list", chrom="list", id="data.frame"),
	  function(object, pos, chrom, id, ...){
		  ## elements of list must be a matrix or an array
		  is.matrix <- all(sapply(object, is, class2="matrix"))
		  is.array <- all(sapply(object, is, class2="array"))
		  stopifnot(is.matrix | is.array)
		  resList <- vector("list", length(object))
		  for(i in seq_along(object)){
			  resList <- segment2(object[[i]], pos=pos[[i]],
					      chrom=chrom[[i]], id=id, ...)
		  }
		  res <- stack(RangedDataList(resList))
		  j <- match("sample", colnames(res))
		  if(length(j) == 1) res <- res[, -j]
		  return(res)
	  })



segmentMatrix <- function(object, pos, chrom, verbose=TRUE, DNAcopy.verbose=0, ...){
	stopifnot(class(object) == "matrix")
	fns <- rownames(object)
	if(any(duplicated(pos))){
		marker.index <- seq_len(nrow(object))[!duplicated(pos)]
	} else marker.index <- seq_len(nrow(object))
	pos <- pos[marker.index]
	chrom <- chrom[marker.index]
	arm <- splitByDistance(pos, thr=75e3)
	index.list <- split(seq_along(marker.index), arm)
	iMax <- sapply(split(marker.index, arm), max)
	pMax <- pos[iMax]
	if(verbose){
		##message("\t\tSmoothing outliers and running circular binary segmentation on ", ncol(object), " samples.")
		pb <- txtProgressBar(min=0, max=length(index.list), style=3)
	}
	md.segs <- vector("list", length(index.list))
	##
	## The names returned by DNAcopy's segment are not necessarily the
	## same as the original names
	## (e.g., if '@' is embedded in the name, or the
	##  first charcter is numeric)
	##  The hash table guarantees that the sample names returned are identical
	##  to the original names
	##
	hash.matrix <- cbind(paste("s", seq_len(ncol(object)), sep=""), colnames(object))
	colnames(hash.matrix) <- c("key", "original.id")
	##
	##
	for(i in seq_along(index.list)){
		if (verbose) setTxtProgressBar(pb, i)
		j <- index.list[[i]]
		CNA.object <- CNA(genomdat=object[j, , drop=FALSE],
				  chrom=chrom[j],
				  maploc=pos[j],
				  data.type="logratio",
				  sampleid=hash.matrix[, "key"])
		smu.object <- smooth.CNA(CNA.object)
		tmp <- segment(smu.object, verbose=DNAcopy.verbose, ...)
		rm(smu.object); gc()
		df <- tmp$output
		sr <- tmp$segRows
		##df <- cbind(tmp$output, tmp$segRows)
		##md.segs[[i]] <-
		firstMarker <- rownames(CNA.object)[sr$startRow]
		endMarker <- rownames(CNA.object)[sr$endRow]
		df$start.index <- match(firstMarker, fns)
		df$end.index <- match(endMarker, fns)
		## if the last marker was duplicated or
		## missing, this might not be true
		stopifnot(max(df$end.index) <= iMax[i])
		md.segs[[i]] <- df
		rm(tmp, df, firstMarker, endMarker, CNA.object); gc()
	}
	if(verbose) close(pb)
	if(length(md.segs) > 1){
		md.segs <- do.call("rbind", md.segs)
	} else md.segs <- md.segs[[1]]
	key.index <- match(md.segs$ID, hash.matrix[, "key"])
	orig.id <- hash.matrix[key.index, "original.id"]
	ranges <- RangedDataCBS(ranges=IRanges(md.segs$loc.start, md.segs$loc.end),
				chromosome=md.segs$chrom,
				sampleId=orig.id,
				coverage=md.segs$num.mark,
				seg.mean=md.segs$seg.mean,
				startIndexInChromosome=md.segs$start.index,
				endIndexInChromosome=md.segs$end.index)
	return(ranges)
}

##.harmonizeSampleNames <- function(original.names, segment.ids){
##	sns <- split(segment.ids, segment.ids)
##	L <- sapply(sns, length)
##	nms <- rep(original.names, L)
##	return(nms)
##}

setMethod("xsegment", signature(object="TrioSet", pedigreeData="Pedigree"),
	  function(object, pedigreeData, id, segment.mindist=TRUE, ..., verbose=FALSE, DNAcopy.verbose=0){
		  ## needs to be ordered
		  if(verbose) message("Segmenting chromosome ", unique(chromosome(object)))
		  ix <- order(chromosome(object), position(object))
		  stopifnot(all(diff(ix) > 0))
		  stopifnot(length(unique(chromosome(object))) == 1)
		  ##
		  ##
		  ##dfl <- vector("list", 22) ## data.frame list
		  if(segment.mindist){
			  if(missing(id)) id <- sampleNames(object)
			  sample.index <- match(id, sampleNames(object))
		  } else{
			  ## all the samples in the samplesheet
			  ##if(missing(id)) id <- allNames(pedigreeData)
			  if(missing(id)){
				  sample.index <- seq_len(nrow(trios(pedigreeData)))
			  } else{
				  trio.index <- trioIndex(pedigreeData)
				  sample.index <- trio.index[match(id, trio.index$individualId), "index.in.pedigree"]
			  }
			  pd <- trios(pedigreeData)[sample.index, ]
			  member.id <- c(t(pd))
		  }
		  stopifnot(length(sample.index) > 0)
		  is.ff <- is(mindist(object), "ff")
		  if(segment.mindist & is.ff){
			  open(mindist(object))
		  } else {
			  if(!segment.mindist & is.ff)
				  open(lrr(object))
		  }
		  fns <- featureNames(object)
		  ##
		  ##
		  ##marker.index.list <- split(seq(length=nrow(minDistanceSet)), chromosome(minDistanceSet))
		  ##for(CHR in 1:22){
		  marker.index <- seq(length=nrow(object))[!duplicated(position(object))]
		  pos <- position(object)[marker.index]
		  chrom <- chromosome(object)[marker.index]
		  if(segment.mindist){
			  CN <- mindist(object)[marker.index, sample.index, drop=FALSE]
			  CN <- matrix(as.numeric(CN), nrow(CN), ncol(CN))
			  dimnames(CN) <- list(featureNames(object)[marker.index], sampleNames(object)[sample.index])
		  } else{
			  ##segment offspring copy number
			  ##offspring.id <- match(id,
			  ##CN <- matrix(NA, length(marker.index), length(sample.index)*3)
			  cnlist <- vector("list", length(sample.index))
			  for(i in seq_along(sample.index)){
				  cnlist[[i]] <- logR(object)[marker.index, sample.index[i], ]
			  }
			  if(length(cnlist) == 1) {
				  CN <- cnlist[[1]]
			  } else CN <- do.call("cbind", cnlist)
			  dimnames(CN) <- list(featureNames(object)[marker.index],
					       member.id)
		  }
		  ##NR <- nrow(object)
		  ##if(verbose) message("Running CBS by chromosome arm")
		  arm <- splitByDistance(pos, thr=75e3)
		  index.list <- split(seq_along(marker.index), arm)
		  iMax <- sapply(split(marker.index, arm), max)
		  pMax <- pos[iMax]
		  md.segs <- list()
		  for(i in seq_along(index.list)){
			  j <- index.list[[i]]
			  CNA.object <- CNA(genomdat=CN[j, , drop=FALSE],
					    chrom=chrom[j],
					    maploc=pos[j],
					    data.type="logratio",
					    sampleid=colnames(CN))
			  smu.object <- smooth.CNA(CNA.object)
			  tmp <- segment(smu.object, verbose=DNAcopy.verbose, ...)
			  df <- tmp$output
			  sr <- tmp$segRows
			  ##df <- cbind(tmp$output, tmp$segRows)
			  ##md.segs[[i]] <-
			  firstMarker <- rownames(CNA.object)[sr$startRow]
			  endMarker <- rownames(CNA.object)[sr$endRow]
			  df$start.index <- match(firstMarker, fns)
			  df$end.index <- match(endMarker, fns)
			  ## if the last marker was duplicated or
			  ## missing, this might not be true
			  stopifnot(max(df$end.index) <= iMax[i])
			  md.segs[[i]] <- df
		  }
		  if(length(md.segs) > 1){
			  md.segs <- do.call("rbind", md.segs)
		  } else md.segs <- md.segs[[1]]
		  if(is.ff){
			  close(mindist(object))
		  }
		  if(!segment.mindist){
			  if(is.ff){
				  close(logR(object))
			  }
		  }
		  ## add 'family' column to ranges
		  ranges <- RangedDataCBS(ranges=IRanges(md.segs$loc.start, md.segs$loc.end),
					  chromosome=md.segs$chrom,
					  sampleId=md.segs$ID,
					  coverage=md.segs$num.mark,
					  seg.mean=md.segs$seg.mean,
					  startIndexInChromosome=md.segs$start.index,
					  endIndexInChromosome=md.segs$end.index)
##		  family <- getFamilyName(ranges, object)
##		  ranges$family <- family
		  if(segment.mindist){
			  mads <- mad.mindist(object)
			  ix <- match(sampleNames(ranges), sampleNames(object))
			  ranges$mindist.mad <- mads[ix]
		  }
		  return(ranges)
})




setMethod("mad.mindist", signature(x="TrioSet"),
	  function(x){
		  x$mindist.mad
	  })

setReplaceMethod("mad.mindist", signature(x="TrioSet"),
		 function(x, value){
			 ## store in phenodata
			 x$mindist.mad <- value
			 return(x)
		 })



##setMethod("xsegment", signature(object="TrioSet"),
##xsegment2 <- function(object, id, ..., verbose=FALSE){
##		  ## needs to be ordered
##		  if(verbose) message("Segmenting chromosome ", unique(chromosome(object)))
##		  ix <- order(chromosome(object), position(object))
##		  stopifnot(all(diff(ix) > 0))
##		  stopifnot(length(unique(chromosome(object))) == 1)
##		  ##
##		  ##
##		  ##dfl <- vector("list", 22) ## data.frame list
####		  if(missing(id)) id <- sampleNames(object)
####		  sample.index <- match(id, sampleNames(object))
####		  stopifnot(length(sample.index) > 0)
####		  open(mindist(object))
##		  fns <- featureNames(object)
##		  ##
##		  ##
##		  ##marker.index.list <- split(seq(length=nrow(minDistanceSet)), chromosome(minDistanceSet))
##		  ##for(CHR in 1:22){
##		  marker.index <- seq(length=nrow(object))[!duplicated(position(object))]
##		  pos <- position(object)[marker.index]
##		  chrom <- chromosome(object)[marker.index]
##		  ##CN <- mindist(object)[marker.index, sample.index, drop=FALSE]
##		  CN <- cbind(logR(object)[marker.index, , 1],
##			      logR(object)[marker.index, , 2],
##			      logR(object)[marker.index, , 3])
##		  colnames(CN) <- c(fatherNames(object),
##				    motherNames(object),
##				    offspringNames(object))
##		  CN <- CN[, !duplicated(colnames(CN))]
##		  ##CN <- matrix(as.numeric(CN), nrow(CN), ncol(CN))
##		  ##dimnames(CN) <- list(featureNames(object)[marker.index], sampleNames(object)[sample.index])
##		  arm <- getChromosomeArm(chrom, pos)
##		  index.list <- split(seq_along(marker.index), arm)
##		  iMax <- sapply(split(marker.index, arm), max)
##		  pMax <- pos[iMax]
##		  cbs.segs <- list()
##		  ##NR <- nrow(object)
##		  ##if(verbose) message("Running CBS by chromosome arm")
##		  for(i in seq_along(index.list)){
##			  j <- index.list[[i]]
##			  CNA.object <- CNA(genomdat=CN[j, , drop=FALSE],
##					    chrom=chrom[j],
##					    maploc=pos[j],
##					    data.type="logratio",
##					    sampleid=colnames(CN))
##			  smu.object <- smooth.CNA(CNA.object)
##			  tmp <- segment(smu.object, verbose=as.integer(verbose), ...)
##			  df <- tmp$output
##			  sr <- tmp$segRows
##			  ##df <- cbind(tmp$output, tmp$segRows)
##			  ##md.segs[[i]] <-
##			  firstMarker <- rownames(CNA.object)[sr$startRow]
##			  endMarker <- rownames(CNA.object)[sr$endRow]
##			  df$start.index <- match(firstMarker, fns)
##			  df$end.index <- match(endMarker, fns)
##			  ## if the last marker was duplicated or
##			  ## missing, this might not be true
##			  stopifnot(max(df$end.index) <= iMax[i])
##			  cbs.segs[[i]] <- df
##		  }
##		  if(length(cbs.segs) > 1){
##			  cbs.segs <- do.call("rbind", cbs.segs)
##		  } else cbs.segs=cbs.segs[[1]]
##		  ##close(mindist(object))
##		  return(cbs.segs)
##}

setMethod("computeBayesFactor", signature(object="TrioSet"),
	  function(object,
		   ranges,
		   returnEmission=FALSE,
		   collapseRanges=TRUE,
		   verbose=TRUE, ...){
		  ## a TrioSet has only one chromosome
		  CHR <- unique(chromosome(object))
		  ranges <- ranges[chromosome(ranges) == CHR, ]
		  stopifnot("pedigreeData" %in% names(list(...)))
		  ##pedigreeData <- list(...)[["pedigreeData"]]
		  if("id" %in% names(list(...))){
			  id <- list(...)[["id"]]
			  ranges <- ranges[sampleNames(ranges) %in% id, ]
		  } else {
			  id <- unique(sampleNames(ranges))
		  }
		  ranges$lik.state <- NA
		  ranges$argmax <- NA
		  ranges$lik.norm <- NA
		  if(verbose){
			  message("\t\tComputing Bayes factors for ", length(id), " files.")
			  pb <- txtProgressBar(min=0, max=length(id), style=3)
		  }
		  for(i in seq_along(id)){
			  if (verbose) setTxtProgressBar(pb, i)
			  this.id <- id[i]
			  k <- match(this.id, sampleNames(object))
			  if(verbose){
				  if(i %% 100 == 0)
					  message("   sample ", this.id, " (", i, "/", length(id), ")")
			  }
			  j <- which(ranges$id == this.id)
			  rd <- joint4(trioSet=object[, k],
				       ranges=ranges[j, ],
				       ...)
			  if(returnEmission) return(rd)
			  ranges$lik.state[j] <- rd$lik.state
			  ranges$argmax[j] <- rd$argmax
			  ranges$lik.norm[j] <- rd$lik.norm
		  }
		  if(verbose) close(pb)
		  if(collapseRanges)
			  ranges <- pruneByFactor(ranges, f=ranges$argmax, verbose=verbose)
		  ranges
	  })

setMethod("todf", signature(object="TrioSet", rangeData="RangedData"),
	  function(object, rangeData, frame, ...){
		  stopifnot(nrow(rangeData) == 1)
		  if(missing(frame)){
			  w <- width(rangeData)
			  frame <- w/0.05  * 1/2
		  }
		  ##overlaps <- findOverlaps(object, rangeData, max.gap=frame)
		  overlaps <- findOverlaps(rangeData, featureData(object), maxgap=frame)
		  marker.index <- matchMatrix(overlaps)[, 2]
		  ##marker.index <- featuresInRangeData(object, rangeData, FRAME=frame)
		  id <- rangeData$id
		  sample.index <- match(id, sampleNames(object))
		  stopifnot(length(sample.index)==1)
		  is.ff <- is(logR(object), "ff")
		  if(is.ff){
			  open(baf(object))
			  open(logR(object))
			  open(mindist(object))
		  }
		  b <- baf(object)[marker.index, sample.index, ]
		  r <- logR(object)[marker.index, sample.index, ]
		  md <- mindist(object)[marker.index, sample.index]
		  if(is.ff){
			  close(baf(object))
			  close(logR(object))
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
	  function(object, ranges, ...){
		  CHR <- unique(chromosome(object))
		  if(verbose) message("Pruning chromosome ", CHR)
		  if(missing(id)) id <- unique(sampleNames(ranges))
		  index <- which(chromosome(ranges) == CHR & sampleNames(ranges) %in% id)
		  ranges <- ranges[index, ]
		  rdList <- vector("list", length(unique(id)))
		  is.ff <- is(mindist(object), "ff")
		  if(is.ff){
			  open(mindist(object))
		  }
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
			  genomdat <- as.numeric(mindist(object)[, k])
			  ## This function currently returns a RangedData object
			  rdList[[j]] <- pruneMD(genomdat,
						 rd,
						 physical.pos=position(object),  ##fD$position,
						 lambda=lambda,
						 MIN.CHANGE=min.change,
						 SCALE.EXP=scale.exp,
						 MIN.COVERAGE=min.coverage)
		  }
		  if(verbose) close(pb)
		  if(is.ff){
			  close(mindist(object))
		  }
		  if(length(rdList) == 1) {
			  rd <- rdList[[1]]
		  } else {
			  rdList <- rdList[!sapply(rdList, is.null)]
			  ##rdList <- lapply(rdList, function(x) as(x, "RangedData"))
			  rdl <- RangedDataList(rdList)
			  rd <- stack(rdl)
			  ##rd <- as(rd, "RangedDataCNV")
			  ix <- match("sample", colnames(rd))
			  if(length(ix) > 0) rd <- rd[, -ix]
		  }
		  ## This will be of class RangedData
		  return(rd)
	 })

## CIDR_Name is not a general label -- need to generalize these functions
##setMethod("offspringNames", signature(object="TrioSet"), function(object){
##	phenoData2(object)[, "CIDR_Name", "O"]
##})
##
##setMethod("fatherNames", signature(object="TrioSet"), function(object){
##	phenoData2(object)[, "CIDR_Name", "F"]
##})
##
##setMethod("motherNames", signature(object="TrioSet"), function(object){
##	phenoData2(object)[, "CIDR_Name", "M"]
##})

setMethod("offspringNames", signature(object="TrioSet"), function(object){
	phenoData2(object)[, "id", "O"]
})

setReplaceMethod("offspringNames", signature(object="TrioSet", value="character"),
		 function(object, value){
			 phenoData2(object)[, "id", "O"] <- value
			 object
		 })

setMethod("fatherNames", signature(object="TrioSet"), function(object){
	phenoData2(object)[, "id", "F"]
})

setReplaceMethod("fatherNames", signature(object="TrioSet", value="character"),
		 function(object, value){
			 phenoData2(object)[, "id", "F"] <- value
			 object
		 })

setMethod("motherNames", signature(object="TrioSet"), function(object){
	phenoData2(object)[, "id", "M"]
})
setReplaceMethod("motherNames", signature(object="TrioSet", value="character"),
		 function(object, value){
			 phenoData2(object)[, "id", "M"] <- value
			 object
		 })

fmoNames <- function(object){
	tmp <- cbind(fatherNames(object), motherNames(object), offspringNames(object))
	colnames(tmp) <- c("F", "M", "O")
	return(tmp)
}


## perhaps make ... additional args to gpar()
##setMethod("xyplotTrioSet", signature(x="formula",
##				     data="TrioSet",
##				     range="RangedDataCNV"),
##xyplotTrioSet <- function(x, data, range, frame=50e3L, ...){
##		  ## for now
##		  ##if(nrow(range) > 1) frame <- 0L
##		  dfList <- vector("list", nrow(range))
##		  for(i in seq_len(nrow(range))){
##			  rm <- findOverlaps(range[i, ], featureData(data), maxgap=frame) ## RangesMatching
##			  mm <- matchMatrix(rm)
##			  mm.df <- data.frame(mm)
##			  mm.df$featureNames <- featureNames(data)[mm.df$subject]
##			  marker.index <- mm.df$subject
##			  sample.index <- match(sampleNames(range)[i], sampleNames(data))
##			  if(any(is.na(sample.index))) stop("sampleNames in RangedData do not match sampleNames in ", class(data), " object")
##			  sample.index <- unique(sample.index)
##			  data2 <- data[marker.index, sample.index]
##			  mm.df$subject <- match(mm.df$featureNames, featureNames(data2))
##			  df <- as(data2, "data.frame")
##			  df$range <- rep(i, nrow(df))##mm.df$query
##			  dfList[[i]] <- df
##		  }
##		  if(length(dfList) == 1) {
##			  df <- dfList[[1]]
##		  } else{
##			  df <- do.call("rbind", dfList)
##		  }
##		  df$range <- factor(df$range, ordered=TRUE, levels=unique(df$range))
##		  if("return.data.frame" %in% names(list(...))){
##			  return.df <- list(...)[["return.data.frame"]]
##			  if(return.df) return(df)
##		  }
##		  list.x <- as.character(x)
##		  ##
##		  ## differs from xyplot2 in VanillaICE here
##		  ##
##		  if("baf" == list.x[[2]]){
##			  df <- df[df$id != "min dist", ]
##		  }
##		  ##
##		  ## end difference
##		  ##
##		  i <- grep("|", list.x, fixed=TRUE)
##		  if(length(i) > 0){
##			  zvar <- list.x[[i]]
##			  zvar <- strsplit(zvar, " | ", fixed=T)[[1]][[2]]
##			  if(zvar == "range"){
##				  tmp <- tryCatch(df$range <- mm.df$query, error=function(e) NULL)
##			  }
##		  }
##		  if("gt" %in% colnames(df)){
##			  xyplot(x, df,
##				 range=range,
##				 id=id,
##				 gt=df$gt,
##				 is.snp=df$is.snp,
##				 ...)
##		  } else {
##			  xyplot(x, df,
##				 id=id,
##				 range=range,
##				 is.snp=df$is.snp,
##				 ...)
##		  }
##	  }

setMethod("xyplot", signature(x="formula", data="TrioSet"),
	  function(x, data, ...){
		  if("range" %in% names(list(...))){
			  ##xyplotTrioSet(x, data, ...)
			  res <- xyplot2(x, data, ...)
		  } else {
			  callNextMethod()
		  }
	  })

setMethod("trioplot", signature(formula="formula", object="TrioSet", range="RangedDataCNV"),
	  function(formula, object, range, ...){
		  xyplot2(x=formula, data=object, range=range, ...)
	  })


setMethod("phenoData2", signature(object="TrioSet"),
	  function(object) object@phenoData2)

setAs("TrioSet", "data.frame",
      function(from, to){
	      ##cn <- copyNumber(from)
	      stopifnot(ncol(from) == 1)
	      cn <- lrr(from)[, 1, ]
	      md <- as.numeric(mindist(from))
	      ##sns <- paste(sampleNames(from), c("F", "M", "O"), sep="_")
	      sns <- phenoData2(from)[, "sampleNames", ]
	      sns <- matrix(sns, nrow(cn), length(sns), byrow=TRUE)
	      sns <- as.character(sns)
	      ##gt <- calls(from)
	      cn <- as.numeric(cn)
	      is.lrr <- c(rep(1L, length(cn)), rep(0L, length(md)))
	      cn <- c(cn, md)
	      sns <- c(sns, rep("min dist", length(md)))
	      ##gt <- as.integer(gt)
	      bf <- as.numeric(baf(from)[, 1, ])
	      bf <- c(bf, rep(NA, length(md)))
	      ##baf.present <- "baf" %in% ls(assayData(from))
	      gt.present <- "call" %in% ls(assayData(from))
	      if(gt.present){
		      gt <- as.numeric(assayDataElement(from, "call"))
		      gt <- c(gt, rep(NA, length(md)))
	      }
	      x <- rep(position(from)/1e6, ncol(from))
	      x <- c(x, position(from)/1e6)
	      ##x <- rep(position(object)[marker.index], 4)/1e6
	      is.snp <- rep(isSnp(from), ncol(from))
	      is.snp <- c(is.snp, isSnp(from))
	      ##id <- rep(sampleNames(from), each=nrow(from))
	      if(!gt.present){
		      df <- data.frame(x=x, lrr=cn, baf=bf, id=sns,
				       is.snp=is.snp,
				       stringsAsFactors=FALSE,
				       is.lrr=is.lrr)
	      } else {
		      df <- data.frame(x=x, lrr=cn,
				       gt=gt,
				       baf=bf,
				       id=sns,
				       is.snp=is.snp,
				       stringsAsFactors=FALSE,
				       is.lrr=is.lrr)
	      }
	      return(df)
      })







setMethod("order", "TrioSet",
	  function(..., na.last=TRUE, decreasing=FALSE){
		  chromosomePositionOrder(...)
	  })
