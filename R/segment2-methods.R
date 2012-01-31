setMethod("segment2", signature(object="TrioSetList"),
	  function(object, md=NULL, segmentParents=TRUE, verbose=TRUE, ...){
		  segmentTrioSetList(object, md, segmentParents=segmentParents, verbose=verbose, ...)
	  })

setMethod("segment2", signature(object="TrioSet"),
	  function(object, md=NULL, segmentParents=TRUE, verbose=TRUE, ...){
		  segmentTrioSet(object, md=md, segmentParents=segmentParents, verbose=verbose, ...)
	  })

setMethod("segment2", signature(object="list"),
	  function(object, pos, chrom, id=NULL, featureNames, segmentParents=TRUE, verbose=TRUE, ...){
		  ## elements of list must be a matrix or an array
		  segmentList(object, pos, chrom, id, featureNames, segmentParents=segmentParents, verbose=verbose, ...)
	  })

setMethod("segment2", signature(object="matrix"),
	  function(object, pos, chrom, id, featureNames, ...){
		  stopifnot(is(id, "character"))
		  segmentMatrix(object, pos, chrom, id, featureNames, ...)
	  })

setMethod("segment2", signature(object="ff_matrix"),
	  function(object, pos, chrom, id, featureNames, ...){
		  segmentff_matrix(object, pos, chrom, id, featureNames, ...)
		  ##segs <- foreach(i=seq_along(ilist), .packages="MinimumDistance") %dopar% segmentMatrix(object[, ilist[[i]]], pos=pos, chrom=chrom, id=id[ilist[[i]]], featureNames, ...)
	  })

setMethod("segment2", signature(object="arrayORff_array"),
	  function(object, pos, chrom, id, featureNames, segmentParents=TRUE, verbose=TRUE, ...){
		  segmentArray(object, pos, chrom, id, featureNames, segmentParents=segmentParents, verbose=verbose, ...)
	  })



segmentTrioSetList <- function(object, md, segmentParents=TRUE, verbose=TRUE, ...){
	if(is.null(md)){
		if(is.null(getCluster())){
			segs <- foreach(trioset=object,
					.packages="MinimumDistance") %do% {
						segment2(object=trioset,
							 segmentParents=segmentParents,
							 verbose=verbose)
					}
			segs <- stackRangedDataList(segs)
		} else {
			segs <- foreach(trioset=object,
					.packages="MinimumDistance") %dopar% {
						segment2(object=trioset,
							 segmentParents=segmentParents,
							 verbose=verbose)
					   }
			segs <- stackRangedDataList(segs)
		}
	} else {
		stopifnot(length(md) == length(chromosome(object)))
		if(is.null(getCluster())){
			segs <- foreach(trioset=object,
					mdElement=md,
					.packages="MinimumDistance",
					.inorder=FALSE) %do% {
						segment2(object=trioset,
							 md=mdElement,
							 verbose=verbose)
					}
			segs <- stackRangedDataList(segs)
		} else {
			segs <- foreach(trioset=object,
					mdElement=md,
					.packages="MinimumDistance",
					.inorder=FALSE) %dopar% {
						segment2(object=trioset,
							 md=mdElement,
							 verbose=verbose)
					}
			segs <- stackRangedDataList(segs)
		}
	}
	return(segs)
}



segmentTrioSet <- function(object, md, segmentParents, verbose, ...){
	if(is.null(md)){
		segmentArray(object=lrr(object),
			     pos=position(object),
			     chrom=chromosome(object),
			     id=trios(object),
			     featureNames=featureNames(object),
			     segmentParents=segmentParents,
			     verbose=verbose)
	} else {
		if(is(md, "logical")){
			md <- mindist(object)
			stopifnot(!is.null(md))
		}
		segment2(object=md,
			 pos=position(object),
			 chrom=chromosome(object),
			 id=sampleNames(object),
			 featureNames=featureNames(object),
			 verbose=verbose)
	}
}



segmentList <- function(object, pos, chrom, id, featureNames, segmentParents=TRUE, verbose=TRUE, ...){
	##warning("segmentList not well tested!")
	dims <- dim(object[[1]])
	if(length(dims) != 2 && length(dims) != 3)
		stop("Elements of list must be a matrix or an array")
	is.matrix <- ifelse(length(dims) == 2, TRUE, FALSE)
	resList <- vector("list", length(object))
	if(!is.matrix & missing(id)) stop("When elements of the list are arrays, a data.frame for 'id' must be provided")
	if(isPackageLoaded("ff")){
		res <- foreach(obj=object, position=pos, chromosome=chrom, fns=featureNames, .inorder=FALSE, .combine=stackRangedDataList, .packages=c("ff", "MinimumDistance")) %dopar% segment2(object=obj, pos=position, chrom=chromosome, id=id, featureNames=fns, verbose=verbose, ...)
	}  else {
		res <- foreach(obj=object, position=pos, chromosome=chrom, fns=featureNames, .inorder=FALSE, .combine=stackRangedDataList, .packages="MinimumDistance") %dopar% segment2(object=obj, pos=position, chrom=chromosome, id=id, featureNames=fns, verbose=verbose, ...)
	}
	j <- match("sample", colnames(res))
	if(!is.na(j)) res <- res[, -j]
	return(res)
}



segmentff_matrix <- function(object, pos, chrom, id, featureNames, ...){
	open(object)
	##ilist <- splitIndicesByLength(seq_len(ncol(object)), 100)
	nc <- ncol(object)
	segs <- vector("list", nc)
	for(i in seq_len(nc)){
		segs[[i]] <- segmentMatrix(object[, i, drop=FALSE], pos=pos, chrom=chrom, id=id[i], featureNames, ...)
	}
	segs <- stack(RangedDataList(segs))
	return(segs)
}
## Read about nesting foreach
##setMethod("segment2", signature(object="arrayORff_array"),



segmentArray <- function(object, pos, chrom, id, featureNames, segmentParents, verbose, ...){
	##open(object)
	## for ff_arrays need to be careful not to pull in too large of a matrix
	##  -- for now, do sample by sample and parallize different chromosomes
	stopifnot(is(id, "data.frame"))
	nc <- ncol(object)
	segs.o <- segs.f <- segs.m <- vector("list", nc)
	if(segmentParents){
		if(verbose) message("segmenting log R ratios for fathers")
		for(i in seq_len(nc)){
			segs.f[[i]] <- segmentMatrix(as.matrix(object[, i, 1]),
						     pos=pos, chrom=chrom, id=id[i, 1], featureNames=featureNames,
						     verbose=verbose, ...)
		}
		if(verbose) message("segmenting log R ratios for mothers")
		for(i in seq_len(nc)){
			segs.m[[i]] <- segmentMatrix(as.matrix(object[, i, 2]),
						     pos=pos, chrom=chrom, id=id[i, 2], featureNames=featureNames,
						     verbose=verbose, ...)
		}
	}
	if(verbose) message("segmenting log R ratios for offspring")
	for(i in seq_len(nc)){
		segs.o[[i]] <- segmentMatrix(as.matrix(object[, i, 3]),
					     pos=pos,
					     chrom=chrom,
					     id=id[i, 3],
					     featureNames=featureNames,
					     verbose=verbose, ...)
	}
	if(!segmentParents){
		segs <- stack(RangedDataList(segs.o))
	} else {
		segs.f <- stack(RangedDataList(segs.f))
		segs.m <- stack(RangedDataList(segs.m))
		segs.o <- stack(RangedDataList(segs.o))
		segs <- stack(RangedDataList(list(segs.f, segs.m, segs.o)))
	}
	return(segs)
}

##setMethod("segment2", signature(object="array"),
##	  function(object, pos, chrom, id, featureNames, ...){
##		  stopifnot(length(dim(object))==3) ## expects a 3d array
##		  J <- dim(object)[[3]]
##		  resList <- vector("list", J)
##		  if("verbose" %in% names(list(...))){
##			  verbose <- list(...)[["verbose"]]
##		  } else verbose <- FALSE
##		  ##featureNames <- rownames(object)
##		  for(j in seq_len(J)){
##			  if(verbose > 0) message("Processing ", ncol(object), " samples")
##			  obj <- object[, , j, drop=FALSE]
##			  dim(obj) <- c(nrow(object), ncol(object))
##			  ##dimnames(obj) <- list(featureNames, id[, j])
##			  resList[[j]] <- segmentMatrix(obj, pos, chrom, id=id[, j], featureNames, ...)
##			  rm(obj)
##		  }
##		  res <- stack(RangedDataList(resList))
##		  j <- match("sample", colnames(res))
##		  if(length(j) == 1) res <- res[, -j]
##		  return(res)
##	  })




segmentMatrix <- function(object, pos, chrom, id, featureNames, ...){
	stopifnot(class(object) == "matrix")
	##featureNames <- rownames(object)
	if(any(duplicated(pos))){
		marker.index <- seq_len(nrow(object))[!duplicated(pos)]
	} else marker.index <- seq_len(nrow(object))
	pos <- pos[marker.index]
	chrom <- chrom[marker.index]
	arm <- splitByDistance(pos, thr=25e6)
#      	arm <- splitByDistance(pos, thr=75e3)
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
				  chrom=chrom[j],
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
