setMethod("initialize", signature(.Object="TrioSetList"),
	  function(.Object,
		   pedigreeData=Pedigree(),
		   assayDataList=AssayDataList(BAF=BAF, logRRatio=logRRatio),
		   logRRatio=array(NA, dim=c(0,0,3)),
		   BAF=array(NA, dim=dim(logRRatio)),
		   featureDataList=GenomeAnnotatedDataFrameFromList(assayDataList),
		   chromosome=integer(),
		   phenoData=annotatedDataFrameFrom(assayDataList, byrow=FALSE),
		   fatherPhenoData=annotatedDataFrameFrom(assayDataList, byrow=FALSE),
		   motherPhenoData=annotatedDataFrameFrom(assayDataList, byrow=FALSE),
		   ...){
		  callNextMethod(.Object,
				 pedigree=pedigreeData,
				 assayDataList=assayDataList,
				 featureDataList=featureDataList,
				 phenoData=phenoData,
				 fatherPhenoData=fatherPhenoData,
				 motherPhenoData=motherPhenoData,
				 chromosome=chromosome,
				 ...)
	  })

setMethod("updateObject", signature(object="TrioSetList"),
	  function(object, ..., verbose=FALSE){
		  if (verbose) message("updateObject(object = 'TrioSetList')")
		  if(!is(object@featureDataList[[1]], "GenomeAnnotatedDataFrame")){
			  fdlist <- lapply(object@featureDataList, updateObject)
			  object@featureDataList <- fdlist
		  }
		  return(object)
	  })

##setMethod("lapply", signature(X="TrioSetList"),
##	  function(X, FUN, ...){
##		  res <- vector("list", length(X))
##		  for(i in seq_along(X)){
##			  res[[i]] <- FUN(X[[i]], ...)
##		  }
##		  res <- new("TrioSetList",
##			     assayData=
##
##		  return(res)
##	  })


GenomeAnnotatedDataFrameFromList <- function(object, annotationPkg){
	nms <- ls(object)
	elt <- object[[nms[1]]]
	fdlist <- vector("list", length(elt))
	for(i in seq_along(elt)){
		fdlist[[i]] <- GenomeAnnotatedDataFrameFromArray(elt[[i]], annotationPkg)
	}
	return(fdlist)
}



TrioSetList <- function(chromosome=integer(),
			pedigreeData=Pedigree(),
			sample.sheet,
			row.names=NULL,
			lrr, baf,
			featureData,
			cdfname){
	if(!missing(lrr)){
		if(!is(lrr[1,1], "integer")){
			stop("lrr should be a matrix of integers. Use integerMatrix(x, scale=100) for the transformation")
		}
		if(!is(baf[1,1], "integer")){
			stop("baf should be a matrix of integers.  Use integerMatrix(x, scale=1000) for the transformation")
		}
	}
	if(nrow(pedigreeData) > 0 & !(missing(lrr) | missing(baf))){
		if(!missing(sample.sheet)){
			if(is.null(row.names)){
				row.names <- rownames(sample.sheet)
			}
			index <- row.names %in% allNames(pedigreeData)
			sample.sheet <- sample.sheet[index, ]
			row.names <- row.names[index]
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
	if(length(chromosome) > 0){
		if(!all(chromosome %in% 1:22)){
			stop("Only autosomal chromosomes (1, 2, ... , 22) allowed")
		}
		if(any(duplicated(chromosome)))
			stop("duplicated chromosomes present")
	} else {
		if(missing(lrr) & missing(baf))
			return(new("TrioSetList"))
	}
	if(missing(lrr) | missing(baf)){
		lrrlist <- baflist <- lapply(chromosome, function(x) array(NA, dim=c(0,0,3)))
		ad <- AssayDataList(BAF=baflist, logRRatio=lrrlist)
		object <- new("TrioSetList",
			      assayDataList=ad,
			      chromosome=chromosome)
		return(object)
	}
	if(!identical(rownames(lrr), rownames(baf))) stop("rownames of lrr and baf must be identical")
	if(missing(featureData)){
		if(missing(cdfname)) stop("if featureData is not supplied, a valid cdfname must be provided for annotating the markers")
		if(any(is.na(rownames(lrr)))){
			message("Removing rows with NA identifiers from lrr & baf matrices")
			lrr <- lrr[!is.na(rownames(lrr)), ]
			baf <- baf[!is.na(rownames(baf)), ]
		}
		##featureData <- oligoClasses:::featureDataFrom(cdfname)
		featureData <- GenomeAnnotatedDataFrameFrom(lrr, cdfname)
		fD <- featureData[order(chromosome(featureData), position(featureData)), ]
		rm(featureData); gc()
	} else {
		if(!is(featureData, "GenomeAnnotatedDataFrame")) stop("featureData must be a GenomeAnnotatedDataFrame")
		fD <- featureData
	}
	if(length(chromosome) > 0){
		fD <- fD[fD$chromosome%in%chromosome, ]
	}
	if(!is.null(rownames(lrr))){
		is.present <- featureNames(fD) %in% rownames(lrr)
		if(!all(is.present)) fD <- fD[is.present, ]
		index <- match(featureNames(fD), rownames(lrr))
		lrr <- lrr[index, ]
		baf <- baf[index, ]
		if(!all(identical(rownames(lrr), sampleNames(fD))))
			stop("rownames of lrr must be the same as the featureNames for the featureData")
	}
	marker.list <- split(seq_along(sampleNames(fD)), fD$chromosome)
	np <- nrow(trios(pedigreeData))
	trio.names <- array(NA, dim=c(length(offspringNames(pedigreeData)), 1, 3))
	dimnames(trio.names) <- list(offspringNames(pedigreeData), "sampleNames", c("F", "M", "O"))
	trio.names[, "sampleNames", ] <- as.matrix(trios(pedigreeData))
	father.names <- originalNames(fatherNames(pedigreeData))
	mother.names <- originalNames(motherNames(pedigreeData))
	offspring.names <- offspringNames(pedigreeData)
	father.index <- match(father.names, colnames(lrr))
	mother.index <- match(mother.names, colnames(lrr))
	offspring.index <- match(offspring.names, colnames(lrr))
	chromosome <- unique(chromosome(fD))
	fdlist <- baflist <- lrrlist <- vector("list", length(chromosome))
	dns <- list(sampleNames(pedigreeData), c("F", "M", "O"))
	for(i in seq_along(marker.list)){
		## Use the name of the offspring as the name for the trio:
		j <- marker.list[[i]]
		nr <- length(j)
		bafArray <- initializeBigArray("baf", dim=c(nr, np, 3), vmode="integer")
		logRArray <- initializeBigArray("lrr", dim=c(nr, np, 3), vmode="integer")
		dimnames(logRArray)[c(2,3)] <- dimnames(bafArray)[c(2,3)] <- dns
		logRArray[,,"F"] <- lrr[j, father.index]
		logRArray[,,"M"] <- lrr[j, mother.index]
		logRArray[,,"O"] <- lrr[j, offspring.index]
		bafArray[,,"F"] <- baf[j, father.index]
		bafArray[,,"M"] <- baf[j, mother.index]
		bafArray[,,"O"] <- baf[j, offspring.index]
		## For each chromosome, create a TrioSet
		lrrlist[[i]] <- logRArray
		baflist[[i]] <- bafArray
		fdlist[[i]] <- fD[j, ]
	}
	ad <- AssayDataList(logRRatio=lrrlist,
			    BAF=baflist)
	object <- new("TrioSetList", assayDataList=ad,
		      featureDataList=fdlist,
		      chromosome=chromosome,
		      pedigree=pedigreeData,
		      fatherPhenoData=fatherPhenoData,
		      motherPhenoData=motherPhenoData,
		      phenoData=phenoData)
	return(object)
}

TrioSetListLD <- function(path, fnames, ext="", samplesheet, row.names,
			  pedigreeData,
			  featureData,
			  annotationPkg, outdir=ldPath(),
			  ffprefix=""){
	if(!is(pedigreeData, "Pedigree")) stop()
	if(missing(featureData)){
		fD <- GenomeAnnotatedDataFrameFrom(file.path(path, paste(fnames[1], ext, sep="")), annotationPkg)
		fD <- fD[chromosome(fD) < 23 & !is.na(chromosome(fD)), ]
	} else {
		fD <- featureData
		rm(featureData); gc()
	}
	ad <- assayDataListLD(path=path,
			      pedigree=pedigreeData,
			      ext=ext,
			      featureData=fD,
			      ffprefix=ffprefix)
	if(!missing(samplesheet)){
		if(missing(row.names)) stop("if samplesheet is provided, row.names can not be missing.")
		index <- row.names %in% allNames(pedigreeData)
		sample.sheet <- samplesheet[index, ]
		row.names <- row.names[index]
		offsprPhenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE,
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
	} else {
		offsprPhenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE, which="offspring")
		fatherPhenoData <- annotatedDataFrameFrom(pedigreeData, FALSE, which="father")
		motherPhenoData <- annotatedDataFrameFrom(pedigreeData, FALSE, which="mother")
	}
	uchrom <- unique(chromosome(fD))
	uchrom <- uchrom[order(uchrom)]
	featureDataList <- vector("list", length(uchrom))
	for(i in seq_along(uchrom)) {
		tmp <- fD[chromosome(fD) == uchrom[i], ]
		featureDataList[[i]] <- tmp[order(position(tmp)), ]
	}
	object <- new("TrioSetList",
		      assayDataList=ad,
		      featureDataList=featureDataList,
		      chromosome=uchrom,
		      pedigree=pedigreeData,
		      fatherPhenoData=fatherPhenoData,
		      motherPhenoData=motherPhenoData,
		      phenoData=offsprPhenoData)
	return(object)
}


setMethod("featureNames", signature(object="TrioSetList"),
	  function(object){
		  lapply(featureDataList(object), sampleNames)
	  })

setMethod("position", signature(object="TrioSetList"),
	  function(object){
		  lapply(featureDataList(object), position)
	  })

setMethod("isSnp", signature(object="TrioSetList"),
	  function(object){
		  lapply(featureDataList(object), function(x) isSnp)
	  })

setMethod("allNames", signature(object="TrioSetList"), function(object) allNames(pedigree(object)))
setMethod("pedigree", signature(object="TrioSetList"), function(object) object@pedigree)
setMethod("trios", signature(object="TrioSetList"), function(object) trios(pedigree(object)))
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
	nchr <- length(chromosome(object))
	ntrios <- ncol(baf(object)[[1]])
	dm <- c(nchr, ntrios)
	names(dm) <- c("chromosomes", "trios")
	return(dm)
})




setMethod("sampleNames", signature(object="TrioSetList"),
	  function(object) offspringNames(object))
##setReplaceMethod("sampleNames", signature(object="TrioSetList", value="character"),
##		 function(object, value){
##			 object <- lapply(object, function(x, value ){
##				 sampleNames(x) <- value
##				 return(x)
##				 }, value=value)
##			 object <- as(object, "TrioSetList")
##			 return(object)
##	 })

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

setMethod("computeBayesFactor", signature(object="TrioSetList"),
	  function(object, ranges,
		   returnEmission=FALSE,
		   collapseRanges=TRUE, ...){
		  computeBayesFactorTrioSetList(object=object,
						ranges=ranges,
						returnEmission=returnEmission,
						collapseRanges=collapseRanges,
						...)
	  })

computeBayesFactorTrioSetList <- function(object,
					  ranges,
					  returnEmission=FALSE,
					  collapseRanges=TRUE,
					  outdir=ldPath(),
					  ...){
	index <- split(seq_len(nrow(ranges)), chromosome(ranges))
	index <- index[names(index) %in% chromosome(object)]
	object <- object[chromosome(object) %in% names(index)]
	index <- index[match(chromosome(object), names(index))]
	##stopifnot(identical(as.character(chromosome(object)), names(index)))
	##if(!identical(as.character(chromosome(object)), names(index))){
	##	stop("The supplied ranges are split into a list by chromosome and that the names
	##}
	X <- i <- NULL
	packages <- neededPkgs()
	map.segs <- foreach(X=object,
			    i=index,
			    .inorder=FALSE,
			    .combine=stackRangedDataList,
			    .packages=packages) %dopar% {
				    computeBayesFactor(object=X,
						       ranges=ranges[i, ],
						       pedigreeData=pedigree(object),
						       outdir=outdir,
						       ...)
			    }
	map.segs$state <- trioStateNames()[map.segs$argmax]
	return(map.segs)
}


setMethod("assayData", signature(object="TrioSetList"),
	  function(object) assayDataList(object))
setMethod("storageMode", "TrioSetList", function(object) storageMode(assayData(object)))

setMethod("phenoData", signature(object="TrioSetList"),
	  function(object) object@phenoData)
setMethod("offspringPhenoData", signature(object="TrioSetList"),
	  function(object) phenoData(object))
setMethod("fatherPhenoData", signature(object="TrioSetList"),
	  function(object) object@fatherPhenoData)
setMethod("motherPhenoData", signature(object="TrioSetList"),
	  function(object) object@motherPhenoData)

setReplaceMethod("assayData", signature=signature(object="TrioSetList",
			      value="AssayData"),
                 function(object, value) {
			 object@assayDataList <- value
			 object
                 })

setMethod("[", signature(x="TrioSetList"),
	  function(x, i, j, ..., drop=FALSE){
		  ## using 'i' to subset markers does not really make
		  ## sense
		  ##
		  ## Use i to subset the list. example, x[1] is still a TrioSetList, but is one chromosome
		  ##
		  if(!missing(i) & !missing(j)){
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- ad[[elt]][i]
				  tmp <- lapply(tmp, function(x, j) {
					  x[, j, , drop=FALSE]
				  }, j=j)
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@chromosome <- chromosome(x)[i]
			  x@featureDataList <- featureDataList(x)[i]
			  x@pedigree <- pedigree(x)[j, ]
			  x@phenoData <- phenoData(x)[j, ]
			  x@fatherPhenoData <- fatherPhenoData(x)[j, ]
			  x@motherPhenoData <- motherPhenoData(x)[j, ]
		  }
		  if(!missing(i) & missing(j)){
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- ad[[elt]][i]
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@chromosome <- chromosome(x)[i]
			  x@featureDataList <- featureDataList(x)[i]
		  }
		  if(missing(i) & !missing(j)){
			  ad <- assayDataList(x)
			  nms <- ls(ad)
			  for(k in seq_along(nms)){
				  elt <- nms[k]
				  tmp <- lapply(ad[[elt]], function(x, j) {
					  x[, j, , drop=FALSE]
				  }, j=j)
				  x <- assayDataElementReplace(x, elt, tmp)
			  }
			  x@pedigree <- pedigree(x)[j, ]
			  x@phenoData <- phenoData(x)[j, ]
			  x@fatherPhenoData <- fatherPhenoData(x)[j, ]
			  x@motherPhenoData <- motherPhenoData(x)[j, ]
		  }
		  return(x)
	  })



setMethod("[[", signature(x="TrioSetList"),
	  function(x, i, j, ..., exact=TRUE){
		  if(missing(i)) return(x)
		  if(length(i) == 1){
			  lrrs <- lrr(x)[[i]]
			  bafs <- baf(x)[[i]]
			  fdlist <- featureDataList(x)[[i]]
			  x <- new("TrioSet",
				   logRRatio=lrrs,
				   BAF=bafs,
				   phenoData=phenoData(x),
				   fatherPhenoData=fatherPhenoData(x),
				   motherPhenoData=motherPhenoData(x),
				   pedigree=pedigree(x),
				   featureData=featureDataList(x)[[i]])
		  } else {
			  stop("subscript out of bounds")
		  }
	  })

setMethod("show", signature(object="TrioSetList"),
	  function(object){
		  lo <- length(lrr(object))
		  cat(class(object), " of length ", lo, "\n", sep="")
	  })

setMethod("length", signature(x="TrioSetList"), function(x) length(x@chromosome))

##setMethod("minimumDistance", signature(object="TrioSetList"),
##	  function(object, narrow.threshold=0.1, ...){
##		  mads.lrr.sample <- mad2(lrr(object), byrow=FALSE)
##		  mads.lrr.marker <- mad2(lrr(object), byrow=TRUE)
##		  mad.sample(object) <- mads.lrr.sample
##		  mad.marker(object) <- mads.lrr.marker
##		  md <- calculateMindist(object)
##		  mads.md <- mad2(md, byrow=FALSE)
##		  mad.mindist(object) <- mads.md
##		  ## add the minimumDistance to the container.
##		  mindist(object) <- md
##		  return(object)
##	  })
setMethod("calculateMindist", signature(object="TrioSetList"),
	  function(object){
		  AssayDataList(calculateMindist(lrr(object)))
	  })


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

setMethod("assayDataList", signature(object="TrioSetList"),
	  function(object)  object@assayDataList)

setMethod("featureDataList", signature(object="TrioSetList"),
	  function(object)  object@featureDataList)

setMethod("lrr", signature(object="TrioSetList"),
	  function(object){
		  ##lapply(object, lrr)
		  assayDataList(object)[["logRRatio"]]
	  })

setMethod("baf", signature(object="TrioSetList"),
	  function(object){
		  ##lapply(object, baf)
		  assayDataList(object)[["BAF"]]
	  })

setMethod("chromosome", signature(object="TrioSetList"),
	  function(object, as.list=FALSE, ...){
		  ##lapply(object, chromosome)
		  if(!as.list) object@chromosome else chromosomeList(object)
	  })

setMethod("chromosomeList", signature(object="TrioSetList"),
	  function(object){
		  ##lapply(object, chromosome)
		  lrrs <- lrr(object)
		  chrom <- rep(object@chromosome, sapply(lrrs, nrow))
		  split(chrom, chrom)
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

setMethod("varLabels", signature(object="TrioSetList"),
	  function(object) varLabels(phenoData(object)))

setMethod("pData", signature(object="TrioSetList"),
	  function(object) pData(phenoData(object)))

setMethod("$", signature(x="TrioSetList"),
	  function(x, name){
		  eval(substitute(phenoData(x)$NAME_ARG, list(NAME_ARG=name)))
	  })



