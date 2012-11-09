setMethod("initialize", "TrioSet",
	  function(.Object,
		   assayData=assayDataNew(logRRatio=logRRatio, BAF=BAF, ...),
		   phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
		   fatherPhenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
		   motherPhenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
		   annotation=character(),
		   featureData=GenomeAnnotatedDataFrameFrom(assayData, annotation, genome=genome),
		   experimentData=new("MIAME"),
		   protocolData=phenoData[, integer(0)],
		   logRRatio=array(NA, dim=c(0, 0, 3)),
		   BAF=array(NA, dim=c(0,0,3)),
		   pedigree=Pedigree(),
		   mindist=NULL,
		   genome=c("hg19", "hg18"), ...){
		  .Object@pedigree <- pedigree
		  .Object@fatherPhenoData <- fatherPhenoData
		  .Object@motherPhenoData <- motherPhenoData
		  callNextMethod(.Object,
				 assayData=assayData,
				 phenoData=phenoData,
				 ##fatherPhenoData=fatherPhenoData,
				 motherPhenoData=motherPhenoData,
				 featureData=featureData,
				 experimentData=experimentData,
				 annotation=annotation,
				 protocolData=protocolData,
				 pedigree=pedigree,
				 mindist=mindist,
				 genome=match.arg(genome), ...)
	  })

setValidity("TrioSet", function(object){
	ped <- pedigree(object)
	validObject(ped)
	validObject(featureData(object))
	nms <- ls(assayData(object))
	if(!all(c("BAF", "logRRatio") %in% nms)){
		msg <- "BAF and logRRatio are required elements of the assayData"
		return(msg)
	}
	elt <- nms[[1]]
	elt <- assayData(object)[[elt]]
	if(ncol(elt) > 0){
		sns.ped <- sampleNames(ped)
		if(length(sns.ped) != ncol(elt)){
			return("Number of samples in pedigree slot should be the same as the number of columns in the TrioSet object")
		}
	}
	if(!identical(sampleNames(object), sampleNames(phenoData(object)))){
		return("sampleNames of TrioSetList object must be the same as the sampleNames of the phenoData")
	}
	if(!identical(fatherNames(object), originalNames(sampleNames(fatherPhenoData(object))))){
		return("fatherNames of TrioSetList object must be the same as the sampleNames of the fatherPhenoData")
	}
	if(!identical(motherNames(object), originalNames(sampleNames(motherPhenoData(object))))){
		stop("motherNames of TrioSetList object must be the same as the sampleNames of the motherPhenoData")
	}
	if(!is.null(mindist(object))){
		if(!identical(colnames(mindist(object)), sampleNames(object)))
			stop("colnames of mindist matrix must be same as the sampleNames of the TrioSet object")
	}
})

setMethod("updateObject", signature(object="TrioSet"),
	  function(object, ..., verbose=FALSE){
		  if (verbose) message("updateObject(object = 'TrioSetList')")
		  if(!is(featureData(object), "GenomeAnnotatedDataFrame")){
			  featureData(object) <- updateObject(featureData(object))
		  }
		  return(object)
	  })

## TrioSet() function fails when this method is uncommented??
##setMethod("dims", signature(object="TrioSet"),
##	  function(object){
##		  nr <- nrow(object)
##		  nchr <- 1
##		  ntrios <- ncol(baf(object))
##		  dm <- c(nchr, ntrios, nr)
##		  names(dm) <- c("chromosomes", "trios", "features")
##		  return(dm)
##	  })
setMethod("pedigree", signature(object="TrioSet"), function(object) object@pedigree)
##setMethod("sampleSheet", signature(object="TrioSet"), function(object) object@sampleSheet)
##setReplaceMethod("sampleSheet", signature(object="TrioSet"), function(object) {object@sampleSheet)
setMethod("lrr", "TrioSet", function(object) assayDataElement(object, "logRRatio"))
setReplaceMethod("lrr", c("TrioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "logRRatio", value)
	 })

setMethod("baf", "TrioSet",
	  function(object) {
		  assayDataElement(object, "BAF")
	 })
setReplaceMethod("baf", c("TrioSet", "array"),
		 function(object, value) {
			 assayDataElementReplace(object, "BAF", value)
	 })
setReplaceMethod("baf", c("TrioSet", "ff_array"),
		 function(object, value) {
			 assayDataElementReplace(object, "BAF", value)
	 })

setMethod("fatherPhenoData", signature(object="TrioSet"),
	  function(object) object@fatherPhenoData)
setMethod("motherPhenoData", signature(object="TrioSet"),
	  function(object) object@motherPhenoData)
setMethod("offspringPhenoData", signature(object="TrioSet"),
	  function(object) phenoData(object))

TrioSet <- function(pedigreeData=Pedigree(),
		    sample.sheet,
		    row.names=NULL,
		    lrr,
		    baf,
		    featureData,
		    cdfname,
		    drop=TRUE,
		    mindist=NULL, genome=c("hg19", "hg18")){
	if(missing(lrr) | missing(baf)){
		object <- new("TrioSet",
			      pedigree=pedigreeData)
		return(object)
	} else{
		if(ncol(lrr) > 0 & nrow(pedigreeData)==0)
			stop("pedigreeData has zero rows")
	}
	if(!missing(lrr) & !missing(baf)){
		if(!identical(rownames(lrr), rownames(baf)))
			stop("rownames of lrr and baf are not identical")
		if(!identical(dim(lrr), dim(baf)))
			stop("lrr and baf must have the same dimension")
		if(!(is(lrr[1,1], "integer") & is(baf[1,1], "integer"))){
			stop("rr and baf must be integers. Use integerMatrix(x, scale=100) to transform log R ratios and integerMatrix(x, scale=1000) for B allele frequencies")
		}
	}
	if(missing(featureData)){
		if(missing(cdfname)) stop("If featureData is not supplied, a valid cdfname must be provided for feature annotation")
		featureData <- GenomeAnnotatedDataFrameFrom(lrr, cdfname, genome=match.arg(genome))
		fD <- featureData[order(chromosome(featureData), position(featureData)), ]
		rm(featureData); gc()
	} else {
		if(!is(featureData, "AnnotatedDataFrame")) stop("featureData must be an AnnotatedDataFrame or a GenomeAnnotatedDataFrame")
		fD <- featureData
	}
	is.present <- sampleNames(fD) %in% rownames(lrr)
	if(!all(is.present)) fD <- fD[is.present, ]
	if(!is.null(rownames(lrr))){
		index <- match(sampleNames(fD), rownames(lrr))
		if(length(index) == 0) {
			if(!missing(cdfname)){
				msg <- paste("rownames for log R ratios do not match feature ids with annotation package ", cdfname)
				stop(msg)
			}
		}
		lrr <- lrr[index, ]
		baf <- baf[index, ]
		stopifnot(all(identical(rownames(lrr), sampleNames(fD))))
	}
	np <- nrow(trios(pedigreeData))
	trio.names <- array(NA, dim=c(length(offspringNames(pedigreeData)), 1, 3))
	dimnames(trio.names) <- list(offspringNames(pedigreeData), "sampleNames", c("F", "M", "O"))
	trio.names[, "sampleNames", ] <- as.matrix(trios(pedigreeData))
	father.names <- fatherNames(pedigreeData)
	mother.names <- motherNames(pedigreeData)
	offspring.names <- offspringNames(pedigreeData)
	father.index <- match(father.names,
			      colnames(lrr))
	if(length(father.index) == 0) stop("father ids in pedigree do not match any of the column names of the lrr matrix")
	mother.index <- match(mother.names,
			      colnames(lrr))
	if(length(mother.index) == 0) stop("mother ids in pedigree do not match any of the column names of the lrr matrix")
	offspring.index <- match(offspring.names,
				 colnames(lrr))
	if(length(offspring.index) == 0) stop("offspring ids in pedigree do not match any of the column names of the lrr matrix")
	nr <- nrow(lrr)
	np <- length(offspring.names)
	bafArray <- initializeBigArray("baf", dim=c(nr, np, 3), vmode="integer")
	logRArray <- initializeBigArray("lrr", dim=c(nr, np, 3), vmode="integer")
	dimnames(bafArray)[[3]] <- dimnames(logRArray)[[3]] <- c("F", "M", "O")
	logRArray[,,"F"] <- lrr[, father.index]
	logRArray[,,"M"] <- lrr[, mother.index]
	logRArray[,,"O"] <- lrr[, offspring.index]
	bafArray[,,"F"] <- baf[, father.index]
	bafArray[,,"M"] <- baf[, mother.index]
	bafArray[,,"O"] <- baf[, offspring.index]
	if(!drop){
		dimnames(bafArray)[c(1,2)] <- dimnames(logRArray)[c(1,2)] <- list(sampleNames(fD), colnames(lrr)[offspring.index])
	}
	if(nrow(pedigreeData) > 0){
		if(!missing(sample.sheet)){
			if(is.null(row.names)){
				row.names <- rownames(sample.sheet)
			}
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
	object <- new("TrioSet",
		      BAF=bafArray,
		      logRRatio=logRArray,
		      phenoData=phenoData,
		      fatherPhenoData=fatherPhenoData,
		      motherPhenoData=motherPhenoData,
		      pedigree=pedigreeData,
		      featureData=fD,
		      mindist=mindist,
		      genome=match.arg(genome))
}



setMethod("show", signature(object="TrioSet"),
	  function(object){
              cat(class( object ), " (storageMode: ", storageMode(object), ")\n", sep="")
              adim <- dim(object)
	      cat("assayData:\n")
              cat("  element names:",
		  paste(assayDataElementNames(object), collapse=", "), "\n")
	      cat("  dimension:\n")
	      print(adim)
	      cat("genome:", genomeBuild(object), "\n")
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
	callNextMethod(object, value)
})

setMethod("mindist", "TrioSet", function(object) object@mindist)
setReplaceMethod("mindist", signature(object="TrioSet", value="matrix"),
		 function(object, value){
			 object@mindist <- value
			 return(object)
		 })

##setReplaceMethod("mindist", signature(object="TrioSet", value="ff_matrix"),
##		 function(object, value){
##			 object@mindist <- value
##			 return(object)
##		 })
##setReplaceMethod("mindist", signature(object="TrioSet", value="NULL"),
##		 function(object, value){
##			 object@mindist <- value
##			 return(object)
##		 })


##setReplaceMethod("trioNames", signature(object="TrioSet"),
##		 function(object,value){
##			 object <- callNextMethod(object, value)
##			 row.names(object@phenoData2) <- value
##			 object
##		 })

setMethod("dim", "TrioSet", function(x) {
	adim <- Biobase:::assayDataDim(assayData(x))
	names(adim) <- c("Features", "Trios", "Members")
	adim
})
setMethod("ncol", signature(x="TrioSet"), function(x) dim(x)[[2]])

setMethod("trios", signature(object="TrioSet"),
	  function(object){
		  trios(pedigree(object))
	  })


setMethod("[", "TrioSet", function(x, i, j, ..., drop = FALSE) {
	if (missing(drop))
		drop <- FALSE
	##	if(length(list(...)) > 0){
	##		k <- list(...)[[1]]
	##	} else k <- NULL
	##	if (missing(i) && missing(j) && is.null(k)) {
	if(missing(i) && missing(j)){
		return(x)
		##if (length(list(...))!=0)
		##	stop("specify features, trios, or samples to subset; use '",
		##	     substitute(x), "$", names(list(...))[[1]],
		##	     "' to access phenoData variables")
		##return(x)
	}
	if (!missing(j) & missing(i)) {
		phenoData(x) <- phenoData(x)[j,, ..., drop = drop]
		protocolData(x) <- protocolData(x)[j,, ..., drop = drop]
		##x@sampleSheet <- x@sampleSheet[j, , , drop=drop]
		##tmp <- pedigree(x)[j, , drop=drop]
		x@pedigree <- pedigree(x)[j, , drop=drop]
		b <- baf(x)[, j, , drop=drop]
		r <- lrr(x)[, j, , drop=drop]
		x <- assayDataElementReplace(x, "logRRatio", r)
		x <- assayDataElementReplace(x, "BAF", b)
		x@fatherPhenoData <- fatherPhenoData(x)[j, ]
		x@motherPhenoData <- motherPhenoData(x)[j, ]
		if(!is.null(mindist(x))){
			mindist(x) <- mindist(x)[, j, drop=FALSE]
		}
		##mad.sample(x) <- mad.sample(x)[j,,...,drop=drop]
	}
	if (!missing(i) & !missing(j)){
		phenoData(x) <- phenoData(x)[j, ..., drop=drop]
		protocolData(x) <- protocolData(x)[j, ..., drop=drop]
		featureData(x) <- featureData(x)[i, ,..., drop=drop]
		b <- baf(x)[i, j, , drop=drop]
		r <- lrr(x)[i, j, , drop=drop]
		##x@sampleSheet <- x@sampleSheet[j, , , drop=drop]
		x@pedigree <- x@pedigree[j, , drop=drop]
		x <- assayDataElementReplace(x, "logRRatio", r)
		x <- assayDataElementReplace(x, "BAF", b)
		x@fatherPhenoData <- fatherPhenoData(x)[j, ]
		x@motherPhenoData <- motherPhenoData(x)[j, ]
		if(!is.null(mindist(x))){
			mindist(x) <- mindist(x)[i, j, drop=FALSE]
		}
	}
	if(!missing(i) & missing(j)){
		featureData(x) <- featureData(x)[i, ,..., drop=drop]
		b <- baf(x)[i, , , drop=drop]
		r <- lrr(x)[i, , , drop=drop]
		x <- assayDataElementReplace(x, "logRRatio", r)
		x <- assayDataElementReplace(x, "BAF", b)
		if(!is.null(mindist(x))){
			mindist(x) <- mindist(x)[i, , drop=FALSE]
		}
	}
	return(x)
})

setMethod("checkOrder", signature(object="TrioSet"),
	  function(object, verbose=FALSE){
		  oligoClasses:::.checkOrder(object, verbose)
	  })

computeBayesFactorTrioSet <- function(object,
				      ranges,
				      returnEmission=FALSE,
				      collapseRanges=TRUE,
				      outdir=ldPath(), ...){
	## a TrioSet has only one chromosome
	ldPath(outdir)
 	CHR <- unique(chromosome(object))
	if(is(ranges, "RangedDataCNV"))
		ranges <- as(ranges, "GRanges")
	ranges <- ranges[chromosome(ranges) == paste("chr",CHR,sep=""), ]
	if(!any(chromosome(ranges) %in% paste("chr", chromosome(object)[1], sep="")))
		stop("chromosome of TrioSet object is not present in the ranges object")
	elementMetadata(ranges)$lik.state <- NA
	elementMetadata(ranges)$argmax <- NA
	elementMetadata(ranges)$lik.norm <- NA
	elementMetadata(ranges)$state <- NA
	if(!"seg.mean" %in% colnames(values(ranges))){
		elementMetadata(ranges)$seg.mean <- NA
	}
	ranges <- ranges[sampleNames(ranges) %in% sampleNames(object), ]
	id <- unique(sampleNames(ranges))
	message("\t\tComputing Bayes factors for ", length(id), " files.")
	pb <- txtProgressBar(min=0, max=length(id), style=3)
	ntrios <- nrow(pedigree(object))
	for(i in seq_along(id)){
		setTxtProgressBar(pb, i)
		this.id <- id[i]
		k <- match(this.id, sampleNames(object))
		if(i %% 100 == 0)
			message("   sample ", this.id, " (", i, "/", length(id), ")")
		j <- which(sampleNames(ranges) == this.id)
		rd <- joint4(id=this.id,
			     trioSet=object[, k],
			     ranges=ranges[j, ],
			     ntrios=ntrios, ...)
			     ##ntrios=ntrios, returnEmission=TRUE, ...)
		if(returnEmission) return(rd)
		if(length(id) > 1){
			elementMetadata(ranges)$lik.state[j] <- values(rd)$lik.state
			elementMetadata(ranges)$argmax[j] <- values(rd)$argmax
			elementMetadata(ranges)$lik.norm[j] <- values(rd)$lik.norm
			elementMetadata(ranges)$state[j] <- trioStateNames()[values(rd)$argmax]
		} else ranges <- rd
	}
	close(pb)
	if(collapseRanges)
		ranges <- pruneByFactor(ranges, f=values(ranges)$argmax, verbose=FALSE)
	ranges
}

setMethod("computeBayesFactor", signature(object="TrioSet"),
	  function(object,
		   ranges,
		   returnEmission=FALSE,
		   collapseRanges=TRUE, outdir=ldPath(), ...){
		  computeBayesFactorTrioSet(object=object,
					    ranges=ranges,
					    returnEmission=returnEmission,
					    collapseRanges=collapseRanges,
					    outdir=outdir,...)
	  })

setMethod("todf", signature(object="TrioSet", rangeData="RangedData"),
	  function(object, rangeData, frame, ...){
		  ## FIX
		  stop("requires mindist(object)")
		  stopifnot(nrow(rangeData) == 1)
		  if(missing(frame)){
			  w <- width(rangeData)
			  frame <- w/0.05  * 1/2
		  }
		  ##overlaps <- findOverlaps(object, rangeData, max.gap=frame)
		  overlaps <- findOverlaps(rangeData, featureData(object), maxgap=frame)
		  marker.index <- subjectHits(overlaps)
		  ##marker.index <- featuresInRangeData(object, rangeData, FRAME=frame)
		  id <- rangeData$id
		  sample.index <- match(id, sampleNames(object))
		  stopifnot(length(sample.index)==1)
		  is.ff <- is(lrr(object), "ff")
		  if(is.ff){
			  open(baf(object))
			  open(lrr(object))
			  open(mindist(object))
		  }
		  b <- baf(object)[marker.index, sample.index, ]
		  r <- lrr(object)[marker.index, sample.index, ]
		  md <- mindist(object)[marker.index, sample.index]
		  if(is.ff){
			  close(baf(object))
			  close(lrr(object))
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
	  function(object, ranges, md, verbose=TRUE, ...){
		  pruneTrioSet(object=object, ranges=ranges, md=md, verbose=verbose, ...)
	 })

pruneTrioSet <- function(object, ranges, md, verbose=TRUE, ...){
	CHR <- unique(chromosome(object))
	if(verbose) message("Pruning chromosome ", CHR)
	id <- unique(sampleNames(ranges))
	index <- which(chromosome(ranges) == CHR & sampleNames(ranges) %in% id)
	ranges <- ranges[index, ]
	rdList <- vector("list", length(unique(id)))
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
		genomdat <- md[, k]
		##genomdat <- as.numeric(mindist(object)[, k])/100
		## This function currently returns a RangedData object
		rdList[[j]] <- pruneMD(genomdat, rd,  ...)
	}
	if(verbose) close(pb)
	rd <- stack(RangedDataList(rdList))
	rd <- rd[, -ncol(rd)]
	return(rd)
}

##setMethod("offspringNames", signature(object="TrioSet"), function(object){
##	phenoData2(object)[, "id", "O"]
##})
##setReplaceMethod("offspringNames", signature(object="TrioSet", value="character"),
##		 function(object, value){
##			 phenoData2(object)[, "id", "O"] <- value
##			 object
##		 })

setMethod("fatherNames", signature(object="TrioSet"), function(object){
	##phenoData2(object)[, "id", "F"]
	fatherNames(pedigree(object))
})
##setReplaceMethod("fatherNames", signature(object="TrioSet", value="character"),
##		 function(object, value){
##			 phenoData2(object)[, "id", "F"] <- value
##			 object
##		 })
setMethod("motherNames", signature(object="TrioSet"), function(object){
	##phenoData2(object)[, "id", "M"]
	motherNames(pedigree(object))
})
##setReplaceMethod("motherNames", signature(object="TrioSet", value="character"),
##		 function(object, value){
##			 phenoData2(object)[, "id", "M"] <- value
##			 object
##		 })
##fmoNames <- function(object){
##	tmp <- cbind(fatherNames(object), motherNames(object), offspringNames(object))
##	colnames(tmp) <- c("F", "M", "O")
##	return(tmp)
##}

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


##setMethod("phenoData2", signature(object="TrioSet"),
##	  function(object) object@phenoData2)
setMethod("allNames", signature(object="TrioSet"), function(object) allNames(pedigree(object)))

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

setMethod("order", signature(...="TrioSet"),
	  function(..., na.last=TRUE,decreasing=FALSE){
		  x <- list(...)[[1]]
		  chromosomePositionOrder(x)
	  })


setMethod("calculateMindist", signature(object="TrioSet"),
	  function(object, verbose=TRUE, ...){
		  calculateMindist(lrr(object))
	  })

setMethod("gcSubtract", signature(object="TrioSet"),
	  function(object, method=c("speed", "lowess"), trio.index, ...){
		  method <- match.arg(method)
		  gcSubtractTrioSet(object, method=method, trio.index, ...)
	  })

gcSubtractTrioSet <- function(object, method, trio.index, ...){
	if(missing(trio.index)) J <- seq_len(ncol(object)) else J <- trio.index
	if(!"gc" %in% fvarLabels(object)) stop("gc not in fvarLabels")
	if(method=="lowess"){
		for(j in J){
			r <- gcSubtractMatrix(lrr(object)[,j,], gc=fData(object)$gc, pos=position(object), ...)
			lrr(object)[,j,] <- integerMatrix(r, 1)
		}
	} else {
		gcbins <- getGcBin(fData(object)$gc)
		for(j in J){
			r <- gcSubtractSpeed(lrr(object)[, j, ], gcbins=gcbins)
			lrr(object)[, j, ] <- integerMatrix(r, 1)
		}
	}
	object
}

gcSubtractSpeed <- function(r, gcbins){
	r.adj <- r
	nc <- ncol(r)
	for(i in seq_along(gcbins)){
		j <- gcbins[[i]]
		mus <- apply(r[j, , drop=FALSE], 2, mean, na.rm=TRUE)
		mus <- matrix(mus, nrow=length(j), ncol=nc, byrow=TRUE)
		r.adj[j, ] <- r[j, , drop=FALSE] - mus
	}
	return(r.adj)
}

getGcBin <- function(gc){
	cuts <- seq(0, 100, by=1)
	bins <- cut(gc, breaks=cuts)
	gcbins <- split(seq_len(length(gc)), bins)  ## contains indices
	L <- sapply(gcbins, length)
	gcbins <- gcbins[L > 0]
	L <- L[L > 0]
	minL <- 50
	while(any(L < minL)){
		index <- which(L < minL)
		if(any(index == 1)){
			gcbins[[2]] <- c(gcbins[[1]], gcbins[[2]])
			gcbins <- gcbins[-1]
			dropFirst <- TRUE
		} else dropFirst <- FALSE
		if(any(index == length(gcbins))){
			LL <- length(gcbins)
			gcbins[[LL-1]] <- c(gcbins[[LL-1]], gcbins[[LL]])
			gcbins <- gcbins[-LL]
			dropLast <- TRUE
		} else dropLast <- FALSE
		if(!(dropFirst | dropLast)){
			index.mid <- index[index > 1 & index < length(gcbins)]
			gcbins[[min(index.mid)-1]] <- c(gcbins[[min(index.mid)-1]], gcbins[[min(index.mid)]])
			gcbins <- gcbins[-min(index.mid)]
		}
		L <- sapply(gcbins, length)
	}
	return(gcbins)
}
