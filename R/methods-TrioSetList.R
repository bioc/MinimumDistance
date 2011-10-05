##setMethod("lapply", signature(X="TrioSetList"),
##	  function(X, FUN, ...){
##		  x <- as(X, "list")
##		  res <- lapply(x, FUN)
##		  X <- as(X, "TrioSetList")
##	  })

setMethod("initialize", signature(.Object="TrioSetList"),
	  function(.Object,
		   pedigree=new("Pedigree"),
		   sampleSheet=new("SampleSheet")){
		  .Object@pedigree <- pedigree
		  .Object@sampleSheet <- sampleSheet
		  return(.Object)
	  })

setMethod("pedigree", signature(object="TrioSetList"), function(object) object@pedigree)
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

TrioSetList <- function(logR, baf,
			pedigree,
			sampleSheet=new("SampleSheet"),
			featureData,
			chromosome=1:22,
			cdfname){
	##stopifnot(is(trioAnnotation, "TrioAnnotation"))
	##phenoDataArray <- as(trioAnnotation, "array")
	##pedigree <- pedigree(trioAnnotation)
	if(missing(featureData)){
		stopifnot(!missing(cdfname))
		featureData <- oligoClasses:::featureDataFrom(cdfname)
		fD <- featureData[order(featureData$chromosome, featureData$position), ]
	} else {
		stopifnot(is(featureData, "AnnotatedDataFrame"))
		fD <- featureData
	}
	marker.list <- split(sampleNames(fD), fD$chromosome)
	marker.list <- marker.list[1:length(marker.list)%in%chromosome]
	np <- nrow(trios(pedigree))
	trioSetList <- vector("list", length(chromosome))
	names(trioSetList) <- 1:length(chromosome)
	father.index <- match(fatherNames(pedigree),
			      colnames(logR))
	mother.index <- match(motherNames(pedigree),
			      colnames(logR))
	offspring.index <- match(offspringNames(pedigree),
				 colnames(logR))
	.Object <- new("TrioSetList", pedigree=pedigree,
		       sampleSheet=sampleSheet)
	for(i in seq_along(marker.list)){
		## Use the name of the offspring as the name for the trio:
		nr <- length(marker.list[[i]])
		bafArray <- logRArray <- array(NA, dim=c(nr, np, 3))
		dimnames(bafArray) <- dimnames(logRArray) <- list(marker.list[[i]],
								  offspringNames(pedigree),
								  colnames(trios(pedigree)))
		##c("F", "M", "O"))
		logRArray[,,"F"] <- logR[marker.list[[i]], father.index]
		logRArray[,,"M"] <- logR[marker.list[[i]], mother.index]
		logRArray[,,"O"] <- logR[marker.list[[i]], offspring.index]
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
					featureData=fD[index,],
					mindist=NULL,
					annotation=cdfname)
		##trioSetList[[chrom]]@phenoData2 <- phenoDataArray
	}
	names(.Object@.Data) <- names(marker.list)
	##trioSetList <- as(trioSetList, "TrioSetList")
	stopifnot(validObject(.Object))
	return(.Object)
}

##setMethod("names", signature(x="TrioSetList") names(x@.Data))

setMethod("mad", signature(x="TrioSetList"), function(x) mad(x[[1]]))

setReplaceMethod("mad", signature(x="TrioSetList", value="ANY"),
		 function(x, value){
			 for(i in seq_along(x)){
				 mad(x[[i]]) <- value
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

setReplaceMethod("mindist", signature(object="TrioSetList"),
		 function(object, value){
			 for(i in seq_along(object)){
				 mindist(object[[i]]) <- value[[i]]
			 }
			 return(object)
		 })

setMethod("minimumDistanceMad", signature(object="TrioSet"),
	  function(object){
		  object$mindist.mad
	  })

setReplaceMethod("minimumDistanceMad", signature(object="TrioSetList"),
		 function(object, value){
			 for(i in seq_along(object)){
				 minimumDistanceMad(object[[i]]) <- value[[i]]
			 }
			 return(object)
		 })

setReplaceMethod("minimumDistanceMad", signature(object="TrioSet"),
		 function(object, value){
			 ## store in phenodata
			 object$mindist.mad <- value
			 return(object)
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
	  function(object){
		  mdList <- lapply(object, calculateMindist)
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
	  function(object, ranges, id, states, baf.sds, mu.logr,
		   log.pi, tau, normal.index, a,
		   prOutlier=c(0.01, 1e-5),
		   prMosaic=0.01,
		   prob.nonMendelian,
		   verbose,
		   returnEmission){
		  if(missing(id)) id <- unique(ranges$id) else stopifnot(id %in% unique(ranges$id))
		  chromosomes <- sapply(object, function(x) unique(chromosome(x)))
		  ranges <- ranges[chromosome(ranges) %in% chromosomes, ]
		  ranges <- ranges[ranges$id %in% id, ]
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
		  for(i in seq_along(object)){
			  if(verbose)
				  message("\tProcessing chromosome ", i, " of ", length(object))
			  CHR <- unique(chromosome(object[[i]]))
			  j <- which(chromosome(ranges) == CHR)
			  if(length(j) < 1) next()
			  rd <- computeBayesFactor(object[[i]],
						   ranges[j, ],
						   pedigreeData=pedigree(object),
						   id=id,
						   states=states,
						   baf.sds=baf.sds,
						   mu.logr=mu.logr,
						   log.pi=log.pi,
						   tau=tau,
						   normal.index=normal.index,
						   a=a,
						   prOutlier=prOutlier,
						   prMosaic=prMosaic,
						   prob.nonMendelian=prob.nonMendelian,
						   returnEmission=returnEmission,
						   verbose=verbose)
			  if(returnEmission) return(rd)
			  ranges$lik.state[j] <- rd$lik.state
			  ranges$argmax[j] <- rd$argmax
			  ranges$lik.norm[j] <- rd$lik.norm
			  ##ranges$DN[j] <- rd$DN
		  }
		  return(ranges)
	  })

setMethod("[", signature(x="TrioSetList"),
	  function(x, i, j, ..., drop=FALSE){
		  if(!missing(i)){
			  nms <- names(x)[i]
			  xlist <- as(x, "list")
			  xlist <- xlist[i]
			  names(xlist) <- nms
			  x <- as(xlist, "TrioSetList")
		  }
		  return(x)
	  })

setMethod("xyplot", signature(x="formula", data="TrioSetList"),
	  function(x, data, ...){
		  stopifnot("range" %in% names(list(...)))
		  range <- list(...)[["range"]]
		  stopifnot(nrow(range)==1)
		  trioSet <- data[[range$chrom]]
		  xyplot(x, trioSet, ...)
	  })

setMethod("chromosome", signature(object="TrioSetList"),
	  function(object) names(object))

setMethod("show", signature(object="TrioSetList"),
	  function(object){
		  lo <- length(object)
		  cat(class(object), " of length ", lo, "\n", sep="")
		  CHR <- chromosome(object)
		  if(lo == 0L){
			  return(NULL)
		  } else {
			  for(i in seq_along(CHR)){
				  adim <- dim(object[[i]])
				  if (length(adim)>1)
					  cat("Chr ", CHR[i], "assayData:",
					      if (length(adim)>1)
					      paste(adim[[1]], "features,",
						    adim[[2]], "samples") else NULL,
					      "\n")
			  }
			  cat(" element names:",
			      paste(assayDataElementNames(object[[1]]), collapse=", "), "\n")
		  }
		  cat("\nSampleSheet:\n")
		  show(sampleSheet(object))
		  cat("\nPedigree:\n")
		  show(pedigree(object))
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
