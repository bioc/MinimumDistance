setMethod("initialize", signature(.Object="TrioSetList"),
	  function(.Object,
		   pedigreeData=new("Pedigree"),
		   sampleSheet=new("SampleSheet")){
		  .Object@pedigree <- pedigreeData
		  .Object@sampleSheet <- sampleSheet
		  return(.Object)
	  })

TrioSetList <- function(lrr, baf,
			pedigreeData,
			sampleSheet,
			featureAnnotation,
			chromosome=1:22,
			cdfname){
	stopifnot(identical(rownames(lrr), rownames(baf)))
	if(missing(featureAnnotation)){
		stopifnot(!missing(cdfname))
		featureAnnotation <- oligoClasses:::featureDataFrom(cdfname)
		fD <- featureAnnotation[order(featureAnnotation$chromosome, featureAnnotation$position), ]
		rm(featureAnnotation); gc()
	} else {
		stopifnot(is(featureAnnotation, "AnnotatedDataFrame"))
		fD <- featureAnnotation
	}
	fD <- fD[order(fD$chromosome, fD$position), ]
	is.present <- sampleNames(fD) %in% rownames(lrr)
	if(!all(is.present)){
		##warning("Excluding SNP ids in featureData not present in rownames of lrr/baf matrices")
		fD <- fD[is.present, ]
	}
	index <- match(sampleNames(fD), rownames(lrr))
	##index <- match(rownames(lrr), sampleNames(fD))
	lrr <- lrr[index, ]
	baf <- baf[index, ]
	stopifnot(all(identical(rownames(lrr), sampleNames(fD))))
	##fD <- fD[index, ]
	if(!missing(sampleSheet)){
		sampleSheet <- sampleSheet[match(allNames(pedigreeData), sampleNames(sampleSheet)), ]
	} else{
		sampleSheet <- SampleSheet(row.names=allNames(pedigreeData))
	}
	marker.list <- split(sampleNames(fD), fD$chromosome)
	##marker.list <- marker.list[1:length(marker.list)%in%chromosome]
	marker.list <- marker.list[names(marker.list)%in%chromosome]
	np <- nrow(trios(pedigreeData))
	trioSetList <- vector("list", length(chromosome))
	names(trioSetList) <- 1:length(chromosome)

	trio.names <- array(NA, dim=c(length(offspringNames(pedigreeData)), 1, 3))
	dimnames(trio.names) <- list(offspringNames(pedigreeData), "sampleNames", c("F", "M", "O"))
	trio.names[, "sampleNames", ] <- as.matrix(trios(pedigreeData))

	father.names <- fatherNames(pedigreeData)
	mother.names <- motherNames(pedigreeData)
	offspring.names <- offspringNames(pedigreeData)
	father.index <- match(father.names,
			      colnames(lrr))
	mother.index <- match(mother.names,
			      colnames(lrr))
	offspring.index <- match(offspring.names,
				 colnames(lrr))
	.Object <- new("TrioSetList",
		       pedigreeData=pedigreeData,
		       sampleSheet=sampleSheet)
	for(i in seq_along(marker.list)){
		## Use the name of the offspring as the name for the trio:
		nr <- length(marker.list[[i]])
		bafArray <- oligoClasses:::initializeBigArray("baf", dim=c(nr, np, 3), vmode="double")
		logRArray <- oligoClasses:::initializeBigArray("lrr", dim=c(nr, np, 3), vmode="double")
		dimnames(bafArray) <- list(marker.list[[i]],
					  sampleNames(pedigreeData),
					  c("F", "M", "O"))
		dimnames(logRArray) <- dimnames(bafArray)
		logRArray[,,"F"] <- lrr[marker.list[[i]], father.index]
		logRArray[,,"M"] <- lrr[marker.list[[i]], mother.index]
		logRArray[,,"O"] <- lrr[marker.list[[i]], offspring.index]
		bafArray[,,"F"] <- baf[marker.list[[i]], father.index]
		bafArray[,,"M"] <- baf[marker.list[[i]], mother.index]
		bafArray[,,"O"] <- baf[marker.list[[i]], offspring.index]
		## For each chromosome, create a TrioSet
		pD <- annotatedDataFrameFrom(as.matrix(logRArray[, , 1]), byrow=FALSE)
		sampleNames(pD) <- colnames(logRArray)
		index <- match(marker.list[[i]], sampleNames(fD))
		## initialize 'TrioSet'
		.Object[[i]] <- new("TrioSet",
				    logRRatio=logRArray,
				    BAF=bafArray,
				    phenoData=pD,
				    phenoArray=trio.names,
				    featureData=fD[index,],
				    annotation=cdfname)
	}
	names(.Object@.Data) <- names(marker.list)
	stopifnot(validObject(.Object))
	return(.Object)
}

setMethod("allNames", signature(object="TrioSetList"), function(object) allNames(pedigree(object)))
setMethod("pedigree", signature(object="TrioSetList"), function(object) object@pedigree)
setMethod("trios", signature(object="TrioSetList"), function(object) trios(pedigree(object)))
setMethod("sampleSheet", signature(object="TrioSetList"), function(object) object@sampleSheet)
setMethod("sampleNames", signature(object="TrioSetList"),
	  function(object) sampleNames(pedigree(object)))
setMethod("nrow", signature(x="TrioSetList"),
	  function(x){
	  sum(sapply(x, nrow))
  })
setMethod("ncol", signature(x="TrioSetList"),
	  function(x) ncol(x[[1]]))
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

setMethod("dims", signature(object="TrioSetList"), function(object){
	names(object) <- paste("chr ", names(object), sep="")
	res <- sapply(object, dim)
	rownames(res)[3] <- c("F, M, O")
	return(res)
})

TrioSetList2 <- function(){


}




##setMethod("names", signature(x="TrioSetList") names(x@.Data))

##setMethod("mad", signature(x="TrioSetList"), function(x) mad(x[[1]]))


setMethod("mindist", signature(object="TrioSetList"), function(object){
	md <- vector("list", length(object))
	for(i in seq_along(object)){
		md[[i]] <- mindist(object[[i]])
	}
	return(md)
})

setReplaceMethod("mindist", signature(object="TrioSetList", value="list"),
		 function(object, value){
			 for(i in seq_along(object)){
				 mindist(object[[i]]) <- value[[i]]
			 }
			 return(object)
		 })


setMethod("order2", "TrioSetList",
	  function(object, ...){
		  orderTrioSetList(object, ...)
	  })

orderTrioSetList <- function(object){
	for(i in seq_along(object)){
		object[[i]] <- order2(object[[i]])
	}
	return(object)
}




setMethod("calculateMindist", signature(object="TrioSetList"),
	  function(object, ..., verbose=TRUE){
		  mdList <- lapply(object, calculateMindist, verbose=verbose)
		  names(mdList) <- paste("chr", names(object))
		  return(mdList)
	  })


setMethod("sampleNames", signature(object="TrioSetList"),
	  function(object) sampleNames(object[[1]]))
setReplaceMethod("sampleNames", signature(object="TrioSetList", value="character"),
		 function(object, value){
			 object <- lapply(object, function(x, value ){
				 sampleNames(x) <- value
				 return(x)
				 }, value=value)
			 object <- as(object, "TrioSetList")
			 return(object)
	 })

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


setMethod("computeBayesFactor", signature(object="TrioSetList", ranges="RangedDataCNV"),
	  function(object, ranges,
		   returnEmission=FALSE, collapseRanges=TRUE, verbose=TRUE, ...){
		  ##if(missing(id)) id <- unique(ranges$id) else stopifnot(id %in% unique(ranges$id))
		  chromosomes <- sapply(object, function(x) unique(chromosome(x)))
		  ranges <- ranges[chromosome(ranges) %in% chromosomes, ]
		  stopifnot(!is.null(mad.marker(object[[1]])))
		  stopifnot(!is.null(mad.sample(object[[1]])))
##		  ranges <- ranges[ranges$id %in% id, ]
##		  if(!"bayes.factor" %in% colnames(ranges)){
##			  ranges$bayes.factor <- NA
##		  }
		  if(!"lik.state" %in% colnames(ranges)){
			  ranges$lik.state <- NA
		  }
		  if(!"lik.norm" %in% colnames(ranges)){
		  	  ranges$lik.norm <- NA
		  }
		  if(!"argmax" %in% colnames(ranges)){
			  ranges$argmax <- NA
		  }
		  if("id" %in% names(list(...))){
			  nsamples <- length(id)
		  } else nsamples <- ncol(object)
		  if("verbose" %in% names(list(...))){
			  verbose <- list(...)[["verbose"]]
		  } else verbose <- FALSE
		  if("returnEmission" %in% names(list(...))){
			  returnEmission <- list(...)[["returnEmission"]]
		  }  else returnEmission <- FALSE
		  if(verbose){
			  message("\t\tComputing Bayes factors for ", length(object), " chromosomes and ", nsamples, " trios.")
			  pb <- txtProgressBar(min=0, max=length(object), style=3)
		  }
		  for(i in seq_along(object)){
			  if (verbose) setTxtProgressBar(pb, i)
##			  if(verbose)
##				  message("\tProcessing chromosome ", i, " of ", length(object))
			  CHR <- unique(chromosome(object[[i]]))
			  j <- which(chromosome(ranges) == CHR)
			  if(length(j) < 1) next()
			  rd <- computeBayesFactor(object[[i]],
						   ranges[j, ],
						   pedigreeData=pedigree(object),
						   collapseRanges=FALSE,
						   verbose=FALSE, ...)
			  if(returnEmission) return(rd)
			  ranges$lik.state[j] <- rd$lik.state
			  ranges$argmax[j] <- rd$argmax
			  ranges$lik.norm[j] <- rd$lik.norm
			  ##ranges$DN[j] <- rd$DN
		  }
		  if(verbose) close(pb)
		  ranges$state <- trioStateNames()[ranges$argmax]
		  if("collapseRanges" %in% names(list(...))){
			  collapseRanges <- list(...)[["collapseRanges"]]
		  } else collapseRanges <- TRUE
		  if(collapseRanges)
			  ranges <- pruneByFactor(ranges, f=ranges$argmax, verbose=verbose)
##		  ranges <- RangedDataMinimumDistance(ranges=ranges(ranges),
##						      values=values(ranges))
		  return(ranges)
 	  })

setMethod("[", signature(x="TrioSetList"),
	  function(x, i, j, ..., drop=FALSE){
		  if(!missing(i) & missing(j)){
			  x@.Data <- x@.Data[i]
		  }
		  if(!missing(i) & !missing(j)){
			  suppressWarnings(x@.Data <- lapply(x, "[", i=i, j=j))
			  x@pedigree <- pedigree(x)[j, ]
		  }
		  if(missing(i) & !missing(j)){
			  suppressWarnings(x@.Data <- lapply(x, "[", j=j))
			  x@pedigree <- pedigree(x)[j, ]
		  }
		  return(x)
	  })



##setMethod("chromosome", signature(object="TrioSetList"),
##	  function(object) names(object))

setMethod("show", signature(object="TrioSetList"),
	  function(object){
		  lo <- length(object)
		  cat(class(object), " of length ", lo, "\n", sep="")
	  })




setMethod("minimumDistance", signature(object="TrioSetList"),
	  function(object, narrow.threshold=0.1, ...){
		  mads.lrr.sample <- mad2(lrr(object), byrow=FALSE)
		  mads.lrr.marker <- mad2(lrr(object), byrow=TRUE)
		  mad.sample(object) <- mads.lrr.sample
		  mad.marker(object) <- mads.lrr.marker
		  md <- calculateMindist(object)
		  mads.md <- mad2(md, byrow=FALSE)
		  mad.mindist(object) <- mads.md
		  ## add the minimumDistance to the container.
		  mindist(object) <- md
		  return(object)
	  })


setMethod("stack", signature(x="TrioSetList"),
	  function(x, ...){
		  bafList=lapply(x, baf)
		  Rs <- sapply(bafList, nrow)
		  C <- ncol(bafList[[1]])
		  logRR <- bf <- array(NA, dim=c(sum(Rs), C, 3))
		  md <- matrix(NA, sum(Rs), C)
		  chrList <- lapply(x, chromosome)
		  chrom <- unlist(chrList)
		  pos <- unlist(lapply(x, position))
		  is.snp <- unlist(lapply(x, isSnp))
		  index <- split(seq_len(sum(Rs)), chrom)
		  for(i in seq_along(x)){
			  j <- index[[i]]
			  bf[j, , ] <- baf(x[[i]])[,,]
			  logRR[j, , ] <- lrr(x[[i]])[,,]
			  md[j, ] <- mindist(x[[i]])[,]
			  ##md.mad[j, ] <- mad(x[[i]])[,]
		  }
		  fns <- as.character(unlist(lapply(x, featureNames)))
		  ##lrr.mad <- apply(md, 2, mad, na.rm=TRUE)
		  dimnames(bf) <- dimnames(logRR) <- list(fns,
							  sampleNames(x[[1]]),
							  c("F","M","O"))
		  featureData <- annotatedDataFrameFrom(as.matrix(bf[,,1]),
							byrow=TRUE)
		  dimnames(md) <- list(fns, sampleNames(x[[1]]))
		  obj <- new("TrioSet",
			     BAF=bf,
			     logRRatio=logRR,
			     mindist=md,
			     phenoData=phenoData(x[[1]]),
			     phenoArray=phenoData2(x[[1]]),
			     featureData=featureData)
		  fData(obj)$chromosome <- chrom
		  fData(obj)$position <- pos
		  fData(obj)$isSnp <- is.snp
		  annotation(obj) <- annotation(x[[1]])
		  return(obj)
	  })

setMethod("lrr", signature(object="TrioSetList"),
	  function(object){
		  lapply(object, lrr)
	  })

setMethod("baf", signature(object="TrioSetList"),
	  function(object){
		  lapply(object, baf)
	  })

setMethod("chromosome", signature(object="TrioSetList"),
	  function(object){
		  lapply(object, chromosome)
	  })
setMethod("position", signature(object="TrioSetList"),
	  function(object){
		  lapply(object, position)
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

