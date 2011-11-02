##setMethod("lapply", signature(X="TrioSetList"),
##	  function(X, FUN, ...){
##		  x <- as(X, "list")
##		  res <- lapply(x, FUN)
##		  X <- as(X, "TrioSetList")
##	  })

setMethod("initialize", signature(.Object="TrioSetList"),
	  function(.Object,
		   pedigreeData=new("Pedigree"),
		   sampleSheet=new("SampleSheet")){
		  .Object@pedigree <- pedigreeData
		  .Object@sampleSheet <- sampleSheet
		  return(.Object)
	  })

setMethod("allNames", signature(object="TrioSetList"), function(object) allNames(pedigree(object)))
setMethod("pedigree", signature(object="TrioSetList"), function(object) object@pedigree)
setMethod("trios", signature(object="TrioSetList"), function(object) trios(pedigree(object)))
setMethod("sampleSheet", signature(object="TrioSetList"), function(object) object@sampleSheet)
setMethod("sampleNames", signature(object="TrioSetList"),
	  function(object) sampleNames(pedigree(object)))
setMethod("nrow", signature(x="TrioSetList"),
	  function(x) nrow(pedigree(x)))
setMethod("offspringNames", signature(object="TrioSetList"), function(object){
	offspringNames(pedigree(object))
})
setMethod("fatherNames", signature(object="TrioSetList"), function(object){
	fatherNames(pedigree(object))
})
setMethod("motherNames", signature(object="TrioSetList"), function(object){
	motherNames(pedigree(object))
})

setMethod("dims", signature(object="TrioSetList"), function(object){
	names(object) <- paste("chr ", names(object), sep="")
	res <- sapply(object, dim)
	rownames(res)[3] <- c("F, M, O")
	return(res)
})

TrioSetList <- function(lrr, baf,
			pedigreeData,
			sampleSheet,
			featureData,
			chromosome=1:22,
			cdfname){
	if(missing(featureData)){
		stopifnot(!missing(cdfname))
		featureData <- oligoClasses:::featureDataFrom(cdfname)
		fD <- featureData[order(featureData$chromosome, featureData$position), ]
		rm(featureData); gc()
	} else {
		stopifnot(is(featureData, "AnnotatedDataFrame"))
		fD <- featureData
	}
	index <- match(rownames(lrr), sampleNames(fD))
	if(any(is.na(index))){
		warning("Some rownames of the log R ratio matrix are not in the corresponding featureData object.")
	}
	fD <- fD[index, ]
	sampleSheet <- sampleSheet[match(allNames(pedigreeData), sampleNames(sampleSheet)), ]
	marker.list <- split(sampleNames(fD), fD$chromosome)
	marker.list <- marker.list[1:length(marker.list)%in%chromosome]
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


##setMethod("names", signature(x="TrioSetList") names(x@.Data))

##setMethod("mad", signature(x="TrioSetList"), function(x) mad(x[[1]]))
setMethod("mad", signature(x="TrioSetList"), function(x) mad(x[[1]]))

setReplaceMethod("mad.sample", signature(x="TrioSetList", value="matrix"),
		 function(x, value){
			 for(i in seq_along(x)){
				 mad.sample(x[[i]]) <- value
			 }
			 return(x)
		 })

setReplaceMethod("mad.marker", signature(x="TrioSetList", value="list"),
		 function(x, value){
			 for(i in seq_along(x)){
				 mad.marker(x[[i]]) <- value[[i]]
			 }
			 return(x)
		 })

setMethod("mindist", signature(object="TrioSetList"), function(object){
	md <- vector("list", length(object))
	for(i in seq_along(object)){
		md[[i]] <- mindist(object[[i]])
	}
	return(md)
})

setReplaceMethod("mindist", signature(object="TrioSetList"),
		 function(object, value){
			 for(i in seq_along(object)){
				 mindist(object[[i]]) <- value[[i]]
			 }
			 return(object)
		 })


setMethod("order", "TrioSetList",
	  function(..., na.last=TRUE, decreasing=FALSE){
		  orderTrioSetList(...)
	  })

orderTrioSetList <- function(object){
	for(i in seq_along(object)){
		object[[i]] <- order(object[[i]])
	}
	return(object)
}





setReplaceMethod("mad.mindist", signature(x="TrioSetList"),
		 function(x, value){
			 for(i in seq_along(x)){
				 mad.mindist(x[[i]]) <- value[[i]]
			 }
			 return(x)
		 })



setMethod("xsegment", signature(object="TrioSetList"),
	  function(object, pedigreeData, id, segment.mindist=TRUE, ...,
		   verbose=FALSE, DNAcopy.verbose=0){
		  ##if(missing(id)) id <- sampleNames(object)
		  ##if(missing(id)) id <- offspringNames(object)
		  dfl <- vector("list", length(object))
		  for(i in seq_along(object)){
			  dfl[[i]] <- xsegment(object[[i]], pedigree(object),
					       id,
					       segment.mindist=segment.mindist,...,
					       verbose=verbose,
					       DNAcopy.verbose=DNAcopy.verbose)
		  }
##		  dfl <- lapply(object, xsegment, pedigreeData=pedigree(object), id,
##				segment.mindist=segment.mindist, ...,
##				verbose=verbose, DNAcopy.verbose=DNAcopy.verbose)
		  ##df <- do.call("rbind", dfl)
		  ranges <- stack(RangedDataList(dfl))
		  index <- match("sample", colnames(ranges))
		  if(length(index) > 0) ranges <- ranges[, -index]
		  return(ranges)
	  })

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
setMethod("ncol", signature(x="TrioSetList"),
	  function(x) ncol(x[[1]]))
setMethod("nrow", signature(x="TrioSetList"),
	  function(x) nrow(x[[1]]))
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



##setMethod("offspringNames", signature(object="TrioSetList"), function(object) offspringNames(object[[1]]))
##setReplaceMethod("offspringNames", signature(object="TrioSetList", value="character"), function(object, value){
##	object <- lapply(object, function(x, value ){
##		offspringNames(x) <- value
##		return(x)
##	}, value=value)
##	object <- as(object, "TrioSetList")
##	return(object)
##})
##setMethod("fatherNames", signature(object="TrioSetList"), function(object) fatherNames(object[[1]]))
##setReplaceMethod("fatherNames", signature(object="TrioSetList", value="character"), function(object, value){
##	object <- lapply(object, function(x, value ){
##		fatherNames(x) <- value
##		return(x)
##	}, value=value)
##	object <- as(object, "TrioSetList")
##	return(object)
##})
##setMethod("motherNames", signature(object="TrioSetList"), function(object) motherNames(object[[1]]))
##setReplaceMethod("motherNames", signature(object="TrioSetList", value="character"), function(object, value){
##	object <- lapply(object, function(x, value ){
##		motherNames(x) <- value
##		return(x)
##	}, value=value)
##	object <- as(object, "TrioSetList")
##	return(object)
##})
##setMethod("fmoNames", signature(object="TrioSetList"), function(object) fmoNames(object[[1]]))

setMethod("computeBayesFactor", signature(object="TrioSetList"),
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
						   returnEmission=returnEmission,
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
		  }
		  if(missing(i) & !missing(j)){
			  suppressWarnings(x@.Data <- lapply(x, "[", j=j))
		  }
		  return(x)
	  })



##setMethod("xyplot", signature(x="formula", data="TrioSetList"),
##	  function(x, data, ...){
##		  stopifnot("rangeData" %in% names(list(...)))
##		  rangeData <- list(...)[["rangeData"]]
##		  stopifnot(nrow(rangeData)==1)
##		  trioSet <- data[[rangeData$chrom]]
##		  xyplotTrioSet(x, data=trioSet, pedigreeData=pedigree(data), ...)
##	  })

setMethod("chromosome", signature(object="TrioSetList"),
	  function(object) names(object))

setMethod("show", signature(object="TrioSetList"),
	  function(object){
		  lo <- length(object)
		  cat(class(object), " of length ", lo, "\n", sep="")
##		  CHR <- chromosome(object)
##		  if(lo == 0L){
##			  return(NULL)
##		  } else {
##			  for(i in seq_along(CHR)){
##				  adim <- dim(object[[i]])
##				  if (length(adim)>1)
##					  cat("Chr ", CHR[i], "assayData:",
##					      if (length(adim)>1)
##					      paste(adim[[1]], "features,",
##						    adim[[2]], "samples") else NULL,
##					      "\n")
##			  }
##			  cat(" element names:",
##			      paste(assayDataElementNames(object[[1]]), collapse=", "), "\n")
##		  }
##		  cat("\nSampleSheet:\n")
##		  show(sampleSheet(object))
##		  cat("\nPedigree:\n")
##		  show(pedigree(object))
	  })

setMethod("show", signature(object="TrioSet"),
	  function(object){
              cat(class( object ), " (storageMode: ", storageMode(object), ")\n", sep="")
              adim <- dim(object)
              if (length(adim)>1)
                  cat("assayData:",
                      if (length(adim)>1)
                      paste(adim[[1]], "features,",
                            adim[[2]], "samples") else NULL,
                      "\n")
              cat("  element names:",
		  paste(assayDataElementNames(object), collapse=", "), "\n")
	  })


setMethod("minimumDistance", signature(object="TrioSetList"),
	  function(object, narrow.threshold=0.1, ...){
		  mads.lrr.sample <- calculateMADlrr(object, by.sample=TRUE)
		  mads.lrr.marker <- calculateMADlrr(object, by.sample=FALSE)
		  mad.sample(object) <- mads.lrr.sample
		  mad.marker(object) <- mads.lrr.marker
		  md <- calculateMindist(object)
		  mads.md <- lapply(md, function(x) apply(x, 2, mad, na.rm=TRUE))
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

setMethod("chromosome", signature(object="TrioSetList"),
	  function(object){
		  lapply(object, chromosome)
	  })
setMethod("position", signature(object="TrioSetList"),
	  function(object){
		  lapply(object, position)
	  })
