TrioSetList <- function(trioAnnotation, logR, baf,
			featureData,
			chromosome=1:22,
			cdfname){
	stopifnot(is(trioAnnotation, "TrioAnnotation"))
	phenoDataArray <- as(trioAnnotation, "array")
	pedigree <- pedigree(trioAnnotation)
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
	np <- nrow(pedigree)
	trioSetList <- vector("list", length(chromosome))
	names(trioSetList) <- 1:length(chromosome)
	father.index <- match(fatherNames(pedigree),
			      colnames(logR))
	mother.index <- match(motherNames(pedigree),
			      colnames(logR))
	offspring.index <- match(offspringNames(pedigree),
				 colnames(logR))
	for(chrom in seq_along(marker.list)){
		## Use the name of the offspring as the name for the trio:
		nr <- length(marker.list[[chrom]])
		bafArray <- logRArray <- array(NA, dim=c(nr, np, 3))
		dimnames(bafArray) <- dimnames(logRArray) <- list(marker.list[[chrom]],
								  offspringNames(pedigree),
								  colnames(pedigree))
		##c("F", "M", "O"))
		logRArray[,,"F"] <- logR[marker.list[[chrom]], father.index]
		logRArray[,,"M"] <- logR[marker.list[[chrom]], mother.index]
		logRArray[,,"O"] <- logR[marker.list[[chrom]], offspring.index]
		bafArray[,,"F"] <- baf[marker.list[[chrom]], father.index]
		bafArray[,,"M"] <- baf[marker.list[[chrom]], mother.index]
		bafArray[,,"O"] <- baf[marker.list[[chrom]], offspring.index]
		## For each chromosome, create a TrioSet
		pD <- annotatedDataFrameFrom(as.matrix(logRArray[, , 1]), byrow=FALSE)
		sampleNames(pD) <- colnames(logRArray)
		index <- match(marker.list[[chrom]], sampleNames(fD))
		## initialize 'TrioSet'
		trioSetList[[chrom]] <- new("TrioSet",
					    logRRatio=logRArray,
					    BAF=bafArray,
					    phenoData=pD,
					    featureData=fD[index,],
					    mindist=NULL,
					    annotation=cdfname)
		## featureData(trioSetList[[chrom]]) <- fD[marker.list[[chrom]], ]
		stopifnot(validObject(trioSetList[[chrom]]))
		trioSetList[[chrom]]@phenoData2 <- phenoDataArray
	}
	trioSetList <- as(trioSetList, "TrioSetList")
	stopifnot(validObject(trioSetList))
	return(trioSetList)
}

setMethod("xsegment", signature(object="TrioSetList"),
	  function(object, id, segment.mindist=TRUE, ..., verbose=FALSE, DNAcopy.verbose=0){
		  if(missing(id)) id <- sampleNames(object)
		  ##if(missing(id)) id <- offspringNames(object)
		  dfl <- lapply(object, xsegment, id=id, segment.mindist=segment.mindist, ..., verbose=verbose,DNAcopy.verbose=DNAcopy.verbose)
		  ##df <- do.call("rbind", dfl)
		  ranges <- stack(RangedDataList(dfl))
		  index <- match("sample", colnames(ranges))
		  if(length(index) > 0) ranges <- ranges[, -index]
		  return(ranges)
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



setMethod("offspringNames", signature(object="TrioSetList"), function(object) offspringNames(object[[1]]))
setReplaceMethod("offspringNames", signature(object="TrioSetList", value="character"), function(object, value){
	object <- lapply(object, function(x, value ){
		offspringNames(x) <- value
		return(x)
	}, value=value)
	object <- as(object, "TrioSetList")
	return(object)
})
setMethod("fatherNames", signature(object="TrioSetList"), function(object) fatherNames(object[[1]]))
setReplaceMethod("fatherNames", signature(object="TrioSetList", value="character"), function(object, value){
	object <- lapply(object, function(x, value ){
		fatherNames(x) <- value
		return(x)
	}, value=value)
	object <- as(object, "TrioSetList")
	return(object)
})
setMethod("motherNames", signature(object="TrioSetList"), function(object) motherNames(object[[1]]))
setReplaceMethod("motherNames", signature(object="TrioSetList", value="character"), function(object, value){
	object <- lapply(object, function(x, value ){
		motherNames(x) <- value
		return(x)
	}, value=value)
	object <- as(object, "TrioSetList")
	return(object)
})
setMethod("fmoNames", signature(object="TrioSetList"), function(object) fmoNames(object[[1]]))

setMethod("computeBayesFactor", signature(object="TrioSetList"),
	  function(object, ranges, id, states, baf.sds, mu.logr,
		   log.pi, tau, normal.index, a,
		   prOutlier=c(0.01, 1e-5),
		   prMosaic=0.01,
		   prob.nonMendelian,
		   df0,
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
						   df0=df0,
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
		  xlist <- as(x, "list")
		  xlist <- xlist[i]
		  x <- as(xlist, "TrioSetList")
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
