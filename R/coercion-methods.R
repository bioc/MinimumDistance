#' Coercion methods in MinimumDistance package
#'
#' @param from see \code{showMethods("coerce")}
#' @param to see \code{showMethods("coerce")}
#' @rdname coercion-methods
#' @aliases coerce,Pedigree,ParentOffspring-method
#' @name coerce
setAs("Pedigree", "ParentOffspring", function(from,to){
  ParentOffspring(id=paste0("trio", nrow(from)),
                  father=fatherNames(from),
                  mother=motherNames(from),
                  offspring=offspringNames(from))
})

#' @aliases coerce,TrioSetList,MinDistExperiment-method
#' @rdname coercion-methods
#' @name coerce
setAs("TrioSetList", "MinDistExperiment", function(from, to){
  trioSet <- stack(from)
  as(trioSet, "MinDistExperiment")
})

#' @aliases coerce,TrioSet,MinDistExperiment-method
#' @rdname coercion-methods
#' @name coerce
setAs("TrioSet", "MinDistExperiment", function(from, to){
  if(ncol(from) > 1) message("only coercing first trio in TrioSet to MinDistExperiment")
  from <- from[, 1]
  ped <- as(pedigree(from), "ParentOffspring")
  gd <- GRanges(paste0("chr", chromosome(from)),
                IRanges(position(from),
                        width=1),
                isSnp=isSnp(from))
  r <- .setColnames(lrr(from)[, 1, ], names(ped))/100
  b <- .setColnames(baf(from)[, 1, ], names(ped))/1000
  assays <- snpArrayAssays(cn=r, baf=b)
  me <- .constructMDE(assays, rowRanges=gd,
                      colData=DataFrame(row.names=names(ped)),
                      ped)
  me
})

#' @aliases coerce,TrioSet,TrioSetList-method
#' @rdname coercion-methods
#' @name coerce
setAs("TrioSet", "TrioSetList",
      function(from, to){
        b <- cbind(baf(from)[, , 1], baf(from)[, , 2], baf(from)[,,3])
        colnames(b) <- c(fatherNames(from),
                         motherNames(from),
                         sampleNames(from))
        r <- cbind(lrr(from)[, , 1], lrr(from)[, , 2], lrr(from)[,,3])
        colnames(r) <- colnames(b)
        TrioSetList(lrr=r,
                    baf=b,
                    pedigreeData=pedigree(from),
                    featureData=featureData(from))
      })

#' @aliases coerce,TrioSet,data.frame-method
#' @rdname coercion-methods
#' @name coerce
setAs("TrioSet", "data.frame",
      function(from, to){
	      ##cn <- copyNumber(from)
	      stopifnot(ncol(from) == 1)
	      cn <- lrr(from)[, 1, ]
	      md <- as.numeric(mindist(from))
	      if(length(md) == 0) stop("minimum distance is not available")
	      ##sns <- paste(sampleNames(from), c("F", "M", "O"), sep="_")
	      ##sns <- phenoData2(from)[, "sampleNames", ]
	      sns <- allNames(from)
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
	      x <- rep(position(from)/1e6, 4)
	      ##x <- c(x, position(from)/1e6)
	      ##x <- rep(position(object)[marker.index], 4)/1e6
	      is.snp <- rep(isSnp(from), 4)
	      ##is.snp <- c(is.snp, isSnp(from))
	      ##id <- rep(sampleNames(from), each=nrow(from))
	      if(!gt.present){
		      df <- data.frame(x=x,
				       lrr=cn,
				       baf=bf,
				       id=sns,
				       is.snp=is.snp,
				       stringsAsFactors=FALSE,
				       is.lrr=is.lrr)
	      } else {
		      df <- data.frame(x=x,
				       lrr=cn,
				       gt=gt,
				       baf=bf,
				       id=sns,
				       is.snp=is.snp,
				       stringsAsFactors=FALSE,
				       is.lrr=is.lrr)
	      }
	      return(df)
            })

#' @param from a \code{TrioSetList}
#' @param to a \code{RangedSummarizedExperiment}
#' @aliases coerce,TrioSetList,RangedSummarizedExperiment-method
#' @rdname TrioSetList-class
setMethod("coerce", signature(from="TrioSetList", to="RangedSummarizedExperiment"),
	  function(from, to){
		  if(ncol(from) > 1) stop("coercion to RangedSummarizedExperiment does not work when ncol > 1")
		  ##nms <- varLabels(from@featureDataList[[1]])
		  chrom <- rep(paste("chr", chromosome(from), sep=""),
			       elementNROWS(from))
		  pos <- unlist(position(from))
		  is.snp <- unlist(lapply(featureDataList(from), isSnp))
		  ## stack the featureDataList to make featureData
		  ## make granges object from featureData
		  sl <- getSequenceLengths(genomeBuild(from))
		  sl <- sl[unique(chrom)]

		  seqinfo <- Seqinfo(seqnames=unique(chrom),
				     genome="hg18")
		  gr <- GRanges(chrom, IRanges(pos,pos), isSnp=is.snp,
				seqlengths=sl,
				seqinfo=seqinfo)
		  names(gr) <- unlist(featureNames(from))
		  rlist <- lrr(from)
		  blist <- baf(from)
		  isff <- is(rlist[[1]], "ff")
		  if(isff) require("ff")
		  ##if(is(rlist[[1]], "ff")
		  rl <- lapply(rlist, "[", , 1, , drop=TRUE) ##function(x) x[, ,drop=FALSE])
		  bl <- lapply(blist, "[", , 1, , drop=TRUE) ##function(x) x[, ,drop=FALSE])
		  r <- do.call("rbind", rl)
		  b <- do.call("rbind", bl)
		  ##rownames(r) <- rownames(b) <- unlist(featureNames(from))
		  ped <- as.character(trios(pedigree(from)))
		  ##colData <- DataFrame(pData(from))
		  ##rownames(colData) <- sampleNames(from)
		  colnames(r) <- colnames(b) <- ped
		  SummarizedExperiment(assays=SimpleList(lrr=r, baf=b),
				       rowRanges=gr)
	  })

#' Coerces a TrioSetList to a TrioSet
#'
#' @param x a \code{TrioSetList}
#' @param ... ignored
#' @return a \code{TrioSet}
#' @aliases stack,TrioSetLiset-method
#' @rdname coercion-methods
#' @export
setMethod("stack", signature(x="TrioSetList"),
	  function(x, ...){
		  b <- baf(x)
		  Rs <- sapply(b, nrow)
		  Cs <- ncol(b[[1]])
		  logRR <- bf <- array(NA, dim=c(sum(Rs), Cs, 3))
		  chrom <- rep(chromosome(x), Rs)
		  ##pos <- unlist(position(x))
		  ##is.snp <- unlist(lapply(x, isSnp))
		  ##is.snp <- unlist(isSnp(x))
		  index <- split(seq_len(sum(Rs)), chrom)
		  for(i in seq_along(x)){
			  j <- index[[i]]
			  bf[j, , ] <- baf(x[[i]])[,,]
			  logRR[j, , ] <- lrr(x[[i]])[,,]
		  }
		  fns <- unlist(featureNames(x))
		  dimnames(bf) <- dimnames(logRR) <- list(fns,
							  sampleNames(x[[1]]),
							  c("F","M","O"))
		  pos <- unlist(position(x))
		  issnp <- unlist(lapply(x@featureDataList, isSnp))
		  featureData <- new("GenomeAnnotatedDataFrame",
				     position=pos,
				     chromosome=chrom,
				     isSnp=issnp,
				     row.names=fns)
		  obj <- new("TrioSet",
			     BAF=bf,
			     logRRatio=logRR,
			     featureData=featureData,
			     pedigree=pedigree(x),
			     motherPhenoData=motherPhenoData(x),
			     fatherPhenoData=fatherPhenoData(x),
			     phenoData=phenoData(x))
		  return(obj)
	  })
