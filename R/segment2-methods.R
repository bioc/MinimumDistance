segmentTrioSetList <- function(object, md, segmentParents=TRUE, verbose=TRUE, ...){
  pkgs <- c("DNAcopy", neededPkgs())
  if(is.null(md)){
    trioset <- NULL
    segs <- foreach(trioset=object,
                    .packages=pkgs, .inorder=FALSE) %dopar% {
                      segment2(object=trioset,
                               segmentParents=segmentParents,
                               verbose=verbose, ...)
                    }
    segs <- stackRangedDataList(segs)
  } else {
    stopifnot(length(md) == length(chromosome(object)))
    trioset <- mdElement <- NULL
    segs <- foreach(trioset=object,
                    mdElement=md,
                    .packages=pkgs,
                    .inorder=FALSE) %dopar% {
                      segment2(object=trioset,
                               md=mdElement,
                               verbose=verbose, ...)
                    }
    segs <- stackRangedDataList(segs)
  }
  metadata(segs) <- list(genome=genomeBuild(object))
  return(segs)
}

segmentTrioSet <- function(object, md, segmentParents, verbose, ...){
  if(is.null(md)){
    segs <- segmentArray(object=lrr(object),
                         pos=position(object),
                         chrom=chromosome(object),
                         id=trios(object),
                         featureNames=featureNames(object),
                         segmentParents=segmentParents,
                         verbose=verbose,
                         genome=genomeBuild(object))
  } else {
    if(is(md, "logical")){
      md <- mindist(object)
      stopifnot(!is.null(md))
    }
    segs <- segment2(object=md,
                     pos=position(object),
                     chrom=chromosome(object),
                     id=sampleNames(object),
                     featureNames=featureNames(object),
                     verbose=verbose,
                     genome=genomeBuild(object))
  }
  metadata(segs) <- list(genome=genomeBuild(object))
  segs
}



segmentList <- function(object, pos, chrom, id, featureNames, segmentParents=TRUE, verbose=TRUE, genome, ...){
  ##warning("segmentList not well tested!")
  dims <- dim(object[[1]])
  if(length(dims) != 2 && length(dims) != 3)
    stop("Elements of list must be a matrix or an array")
  is.matrix <- ifelse(length(dims) == 2, TRUE, FALSE)
  resList <- vector("list", length(object))
  if(!is.matrix & missing(id)) stop("When elements of the list are arrays, a data.frame for 'id' must be provided")
  pkgs <- unique(c("DNAcopy", neededPkgs()))
  if(isPackageLoaded("ff")){
    pkgs <- c("ff", pkgs)
    obj <- position <- chromosome <- fns <- NULL
    res <- foreach(obj=object,
                   position=pos,
                   chromosome=chrom,
                   fns=featureNames, .inorder=FALSE,
                   ##.combine=stackRangedDataList,
                   .packages=pkgs) %dopar% segment2(object=obj, pos=position, chrom=chromosome, id=id, featureNames=fns, verbose=verbose, genome=genome, ...)
  }  else {
    obj <- position <- chromosome <- fns <- NULL
    res <- foreach(obj=object, position=pos, chromosome=chrom, fns=featureNames, .inorder=FALSE, .packages=pkgs) %dopar% segment2(object=obj, pos=position, chrom=chromosome, id=id, featureNames=fns, verbose=verbose, genome=genome, ...)
  }
  res <- unlist(GRangesList(res))
  ##	j <- match("sample", colnames(res))
  ##	if(!is.na(j)) res <- res[, -j]
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
  ##segs <- stack(RangedDataList(segs))
  segs <- unlist(GRangesList(segs))
  return(segs)
}

segmentff_matrix2 <- function(object, pos, chrom, fid, mid, oid,
			      featureNames, sample.index, ...){
  segs.f <- vector("list", length(sample.index))
  segs.m <- vector("list", length(sample.index))
  segs.o <- vector("list", length(sample.index))
  for(i in seq_along(sample.index)){
    j <- sample.index[i]
    obj <- as.matrix(object[, j, 1])
    segs.f[[i]] <- segmentMatrix(obj, pos=pos,
                                 chrom=chrom,
                                 id=fid[i],
                                 featureNames, ...)
    obj <- as.matrix(object[, j, 2])
    segs.m[[i]] <- segmentMatrix(obj, pos=pos,
                                 chrom=chrom,
                                 id=mid[i],
                                 featureNames, ...)
    obj <- as.matrix(object[, j, 3])
    segs.o[[i]] <- segmentMatrix(obj, pos=pos,
                                 chrom=chrom,
                                 id=oid[i],
                                 featureNames, ...)
  }
  segs.f <- unlist(GRangesList(segs.f))
  segs.m <- unlist(GRangesList(segs.m))
  segs.o <- unlist(GRangesList(segs.o))
  segs <- unlist(GRangesList(segs.f, segs.m, segs.o))
  ##segs <- stack(RangedDataList(segs.f, segs.m, segs.o))
  return(segs)
}

segmentArray <- function(object, pos, chrom, id, featureNames, segmentParents, verbose, ...){
  ##open(object)
  ## for ff_arrays need to be careful not to pull in too large of a matrix
  ##  -- for now, do sample by sample and parallize different chromosomes
  if(!is(id, "data.frame")) stop("'id' should be a data.frame")
  nc <- ncol(object)
  segs <- segs.o <- segs.f <- segs.m <- vector("list", nc)
  if(is(segmentParents, "logical")){
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
      ##segs <- stack(RangedDataList(segs.o))
      ##segs <- stackRangedDataList(segs.o)
      segs <- unlist(GRangesList(segs.o))
    } else {
      ##segs.f <- stackRangedDataList(segs.f)
      segs.f <- unlist(GRangesList(segs.f))
      segs.m <- unlist(GRangesList(segs.m))
      segs.o <- unlist(GRangesList(segs.o))
      ##			segs.m <- stackRangedDataList(segs.m)
      ##			segs.o <- stackRangedDataList(segs.o)
      segs <- unlist(GRangesList(segs.f, segs.m, segs.o))
      ##segs <- stackRangedDataList(list(segs.f, segs.m, segs.o))
    }
  } else {
    if(!segmentParents %in% 1:3) stop("segmentParents must be logical, or an integer (1, 2, or 3)")
    for(i in seq_len(nc)){
      k <- segmentParents
      segs[[i]] <- segmentMatrix(as.matrix(object[, i, k]),
                                 pos=pos,
                                 chrom=chrom,
                                 id=id[i, k],
                                 featureNames=featureNames,
                                 verbose=verbose, ...)
    }
    ##segs <- stackRangedDataList(segs)
    segs <- unlist(GRangesList(segs))
  }
  return(segs)
}


segmentMatrix <- function(object, pos, chrom, id, featureNames,
			  genome, gapsize=75e3, ...){
  if(!is(object, "matrix"))
    stop("object must be a matrix")
  ##featureNames <- rownames(object)
  if(any(duplicated(pos))){
    marker.index <- seq_len(nrow(object))[!duplicated(pos)]
  } else marker.index <- seq_len(nrow(object))
  pos <- pos[marker.index]
  chrom <- chrom[marker.index]
  arm <- splitByDistance(pos, thr=gapsize)
  index.list <- split(seq_along(marker.index), arm)
  iMax <- sapply(split(marker.index, arm), max)
  pMax <- pos[iMax]
  hash.matrix <- cbind(paste("s", seq_len(ncol(object)), sep=""), id)
  colnames(hash.matrix) <- c("key", "original.id")
  rownames(object) <- featureNames
  segs <- vector("list", length(index.list))
  for(i in seq_along(index.list)){
    ##if (verbose) setTxtProgressBar(pb, i)
    j <- index.list[[i]]
    r <- object[j, , drop=FALSE]/100
    CNA.object <- CNA(genomdat=r,
                      chrom=chrom[j],
                      maploc=pos[j],
                      data.type="logratio",
                      sampleid=hash.matrix[, "key"])
    smu.object <- smooth.CNA(CNA.object)
    tmp <- segment(smu.object, ...)
    rm(smu.object); gc()
    df <- tmp$output
    sr <- tmp$segRows
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
  ranges <- GRanges(paste("chr", segs$chrom, sep=""),
                    IRanges(segs$loc.start, segs$loc.end),
                    sample=orig.id,
                    numberProbes=segs$num.mark,
                    seg.mean=segs$seg.mean)
  if(!missing(genome)){
    sl <- getSequenceLengths(genome)
    sl <- sl[match(unique(as.character(seqnames(ranges))), names(sl))]
    seqlengths(ranges) <- sl
  }
  return(ranges)
}


.unmake.names <- function(id, original_id){
  mapping <- setNames(original_id, make.names(original_id))
  as.character(mapping[as.character(id)])
}

.dnacopy2granges <- function(x, seq_info, original_id){
  ## The names in x$ID might not be the same as the original because of make.names.
  ids <- .unmake.names(x$ID, original_id)
  g <- GRanges(x$chrom,
               IRanges(x$loc.start, x$loc.end),
               sample=ids,
               numberProbes=x$num.mark,
               seg.mean=x$seg.mean)
  seqlevels(g) <- seqlevels(seq_info)
  seqinfo(g) <- seq_info
  g
}


.smoothAndSegment <- function(x, rowdata, param){
  CNA.object <- CNA(genomdat=x,
                    chrom=as.character(seqnames(rowdata)),
                    maploc=start(rowdata), ##pos[j],
                    data.type="logratio",
                    sampleid=colnames(x))
  smu.object <- smooth.CNA(CNA.object)
  segment(smu.object, alpha=alpha(param),
          min.width=min.width(param),
          undo.splits=undo.splits(param))$output
}

.setFilename <- function(granges, se){
  files <- se$filename
  mdfiles <- setNames(files[offspring(se)], .mindistnames(offspring(se)))
  files <- c(files, mdfiles)
  files[granges$sample]
}

.validNames <- function(g, id){
  if(!all(g$sample %in% id)){
    warning("appears that names were mangled by dnacopy. Attempting to unmangle")
    g$sample <- gsub("\\.", "-", g$sample)
    if(!all(g$sample %in% id)){
      message("unmangling not successful. Use _ instead of symbols like '-' to avoid name mangling")
      return(FALSE)
    }
  }
  TRUE
}


.pedigreeId <- function(object) {
  paste(colData(object)$filename["father"],
        colData(object)$filename["mother"], sep="_")
}

subsetGRangesById <- function(g, id){ g[g$sample %in% id]}



## Narrow the minimum distance segmentation breakpoints
##
## Narrow the minimum distance segmentation interval if breakpoints in
## the offspring occur within the interval.
##
## @param offspring_grl a \code{GRangesList} object of the offspring segmentation
## @param mindist_grl  a \code{GRangesList} object of the minimum distance segmentation
## @param mads a numeric vector containing the minimum distance median absolute deviations
## @param param
## @export
narrow2 <- function(offspring_grl, mindist_grl, mads, param){
  mindist_grl2 <- foreach(md_gr=mindist_grl, offspr_gr=offspring_grl, md.mad=mads) %do%{
    .narrowMinDistGRanges(md_gr=md_gr, offspr_gr=offspr_gr, md.mad=md.mad, param=param)
  }
  setNames(GRangesList(mindist_grl2), names(mindist_grl))
}

.narrowMinDistGRanges <- function(md_gr, offspr_gr, md.mad, param){
  ## threshold is the number of mads
  object <- md_gr
  keep <- segMeanAboveThr(mean=md_gr$seg.mean, mad=md.mad, nmad=nMAD(param))
  md.below.thr <- md_gr[!keep]
  md_gr <- md_gr[keep]
  if(length(md_gr) < 1) return(object)
  offspr_gr <- subsetByOverlaps(offspr_gr, md_gr)
  disj <- disjoin(c(md_gr, offspr_gr))
  disj <- subsetByOverlaps(disj, md_gr)
  hits <- findOverlaps(md_gr, disj)
  ## which minimumdistance intervals are spanned by a disjoint interval
  j <- subjectHits(hits) ##r
  i <- queryHits(hits)   ##s
  ##disj$sample <- names(object)
  disj$sample <- md_gr$sample[1]
  disj$seg.mean <- NA
  disj$seg.mean[j] <- md_gr$seg.mean[i]
  disj$filename <- md_gr$filename[1]
  ##
  ## filtered minimum distance ranges (md < thr) will only be in the
  ## offspring segments
  ##
  mcols(md.below.thr) <- mcols(md.below.thr)[, -grep("numberProbes", colnames(mcols(md.below.thr)))]
  disj <- sort(c(disj, md.below.thr))
  return(disj)
}
