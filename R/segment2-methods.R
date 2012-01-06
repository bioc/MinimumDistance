
setMethod("segment2", signature(object="list", pos="list", chrom="list"),
	  function(object, pos, chrom, id=NULL, featureNames, ...){
		  ## elements of list must be a matrix or an array
		  is.matrix <- is(object[[1]], "matrix") || is(object[[1]], "ff_matrix")
		  is.array <- class(object[[1]]) == "array" || is(object[[1]], "ff_array")
		  stopifnot(is.matrix | is.array)
		  resList <- vector("list", length(object))
		  for(i in seq_along(object)){
			  if(is.array){
				  if(missing(id)) stop("When elements of the list are arrays, a data.frame for 'id' must be provided")
				  resList[[i]] <- segment2(object[[i]], pos=pos[[i]],
							   chrom=chrom[[i]], id=id,
							   featureNames=featureNames[[i]], ...)
			  } else {
##				  if(is.null(colnames(object[[1]]))){
##					  stop("column names of the elements in the list can not be null.")
##				  }
##				  if(!is.null(id))
##					  message("id is ignored as elements of list are matrices. Column names of the matrices will be used to label the segments")
				  resList[[i]] <- segment2(object[[i]], pos=pos[[i]],
							   chrom=chrom[[i]], id=id,
							   featureNames=featureNames[[i]], ...)
			  }
		  }
		  res <- stack(RangedDataList(resList))
		  j <- match("sample", colnames(res))
		  if(length(j) == 1) res <- res[, -j]
		  return(res)
	  })

setMethod("segment2", signature(object="matrix"),
	  function(object, pos, chrom, id, featureNames, ...){
		  res <- segmentMatrix(object, pos, chrom, id, featureNames, ...)
	  })

setMethod("segment2", signature(object="ff_matrix"),
	  function(object, pos, chrom, id, featureNames, ...){
		  open(object)
		  segment2(object[,], pos=pos, chrom=chrom, id=id, featureNames, ...)
	  })

setMethod("segment2", signature(object="array", pos="ANY", chrom="ANY", id="data.frame"),
	  function(object, pos, chrom, id, featureNames, ...){
		  stopifnot(length(dim(object))==3) ## expects a 3d array
		  J <- dim(object)[[3]]
		  resList <- vector("list", J)
		  if("verbose" %in% names(list(...))){
			  verbose <- list(...)[["verbose"]]
		  } else verbose <- FALSE
		  ##featureNames <- rownames(object)
		  for(j in seq_len(J)){
			  if(verbose > 0) message("Processing ", ncol(object), " samples")
			  obj <- object[, , j, drop=FALSE]
			  dim(obj) <- c(nrow(object), ncol(object))
			  ##dimnames(obj) <- list(featureNames, id[, j])
			  resList[[j]] <- segmentMatrix(obj, pos, chrom, id=id[, j], featureNames, ...)
			  rm(obj)
		  }
		  res <- stack(RangedDataList(resList))
		  j <- match("sample", colnames(res))
		  if(length(j) == 1) res <- res[, -j]
		  return(res)
	  })

setMethod("segment2", signature(object="ff_array", id="data.frame"),
	  function(object, pos, chrom, id, featureNames, ...){
		  open(object)
		  segment2(object[,,], pos=pos, chrom=chrom, id=id, featureNames, ...)
	  })


segmentMatrix <- function(object, pos, chrom, id, featureNames, ...){
	stopifnot(class(object) == "matrix")
	##featureNames <- rownames(object)
	if(any(duplicated(pos))){
		marker.index <- seq_len(nrow(object))[!duplicated(pos)]
	} else marker.index <- seq_len(nrow(object))
	pos <- pos[marker.index]
	##chrom <- chrom[marker.index]
	arm <- splitByDistance(pos, thr=75e3)
	index.list <- split(seq_along(marker.index), arm)
	iMax <- sapply(split(marker.index, arm), max)
	pMax <- pos[iMax]
##	if(verbose){
##		##message("\t\tSmoothing outliers and running circular binary segmentation on ", ncol(object), " samples.")
##		pb <- txtProgressBar(min=0, max=length(index.list), style=3)
##	}
	##
	## The names returned by DNAcopy's segment are not necessarily the
	## same as the original names
	## (e.g., if '@' is embedded in the name, or the
	##  first charcter is numeric)
	##  The hash table guarantees that the sample names returned are identical
	##  to the original names
	##
	##hash.matrix <- cbind(paste("s", seq_len(ncol(object)), sep=""), colnames(object))
	hash.matrix <- cbind(paste("s", seq_len(ncol(object)), sep=""), id)
	colnames(hash.matrix) <- c("key", "original.id")
	rownames(object) <- featureNames
	##
	##
	segs <- vector("list", length(index.list))
	for(i in seq_along(index.list)){
		##if (verbose) setTxtProgressBar(pb, i)
		j <- index.list[[i]]
		CNA.object <- CNA(genomdat=object[j, , drop=FALSE],
				  chrom=rep(chrom, length(j)),
				  maploc=pos[j],
				  data.type="logratio",
				  sampleid=hash.matrix[, "key"])
		smu.object <- smooth.CNA(CNA.object)
		tmp <- segment(smu.object, ...)
		rm(smu.object); gc()
		df <- tmp$output
		sr <- tmp$segRows
		##df <- cbind(tmp$output, tmp$segRows)
		##md.segs[[i]] <-
		firstMarker <- rownames(CNA.object)[sr$startRow]
		endMarker <- rownames(CNA.object)[sr$endRow]
		df$start.index <- match(firstMarker, featureNames)
		df$end.index <- match(endMarker, featureNames)
		## if the last marker was duplicated or
		## missing, this might not be true
		stopifnot(max(df$end.index) <= iMax[i])
		segs[[i]] <- df
		rm(tmp, df, firstMarker, endMarker, CNA.object); gc()

	}
	##if(verbose) close(pb)
	if(length(segs) > 1){
		segs <- do.call("rbind", segs)
	} else segs <- segs[[1]]
	key.index <- match(segs$ID, hash.matrix[, "key"])
	orig.id <- hash.matrix[key.index, "original.id"]
	ranges <- RangedDataCBS(ranges=IRanges(segs$loc.start, segs$loc.end),
				chromosome=segs$chrom,
				sampleId=orig.id,
				coverage=segs$num.mark,
				seg.mean=segs$seg.mean,
				startIndexInChromosome=segs$start.index,
				endIndexInChromosome=segs$end.index)
	return(ranges)
}
