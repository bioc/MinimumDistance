AssayDataList <- function(storage.mode = c("lockedEnvironment", "environment", "list"), ...) {
	storage.mode <- match.arg(storage.mode) ## defaults to "lockedEnvironment"
	assayData <- switch(storage.mode,
			    lockedEnvironment =,
			    environment = new.env(parent=emptyenv()),
			    list = list())
	arglist <- list(...)
	for (nm in names(arglist)) assayData[[nm]] <- arglist[[nm]]
	if (storage.mode == "lockedEnvironment") Biobase:::assayDataEnvLock(assayData)
	assayData
}

assayDataListDims <- function(object) {
	nms <- if (assayDataStorageMode(object) == "list") names(object) else ls(object)
	if (length(nms) == 0)
		return(matrix(integer(0), nrow = 2, ncol = 0,
			      dimnames = list(c("Features", "Samples"), character(0))))
	d <- lapply(nms, function(i) lapply(object[[i]], dim)) ##dim(object[[i]]))
	##rownames(d) <- c("Features", "Samples", rep("...", nrow(d)-2))
	names(d) <- nms
	##colnames(d) <- nms
	##d[,order(colnames(d)), drop=FALSE]
	return(d)
}


validAssayDataDims <- function(object){
	msg <- NULL
	d <- assayDataListDims(object)
	firstElement <- d[[1]]
	d <- d[-1]
	res <- sapply(d, function(i) identical(i, firstElement))
	## check that the 3rd dimension is 3
	if(!all(res)){
		msg <- "Assay data elements must have the same dimension"
	}
	if(firstElement[[1]][3] != 3){
		msg <- c(msg, "third dimension of assayData elements must be 3")
	}
	if(is.null(msg)) return(TRUE) else return(msg)
}

assayDataListLD <- function(path="", ext="", pedigree, featureData){
	filenames <- paste(allNames(pedigree), ext, sep="")
	fnames <- file.path(path, filenames)
	stopifnot(all(file.exists(fnames)))
	if(missing(featureData)) stop("featureData can not be missing")
	if(missing(pedigree)) stop("pedigree can not be missing")
	index <- split(seq_len(nrow(featureData)), chromosome(featureData))
	if(isPackageLoaded("ff")){
		message("Initializing ff arrays for BAFs and LRRs (this may take some time) ...")
	} else {
		message("ff package not loaded -- initializing arrays for BAFs and LRRs ...")
	}
	baflist <- lapply(index, function(x) initializeBigArray("baf", dim=c(length(x), nrow(pedigree), 3), vmode="double"))
	lrrlist <- lapply(index, function(x) initializeBigArray("lrr", dim=c(length(x), nrow(pedigree), 3), vmode="double"))
	if(isPackageLoaded("ff"))
		message("Reading ", length(fnames),
			" files and writing ff files to directory ", ldPath())
	baflist <- lapply(baflist, function(x, sampleNames) {
		colnames(x) <- sampleNames
		return(x)
	}, sampleNames=sampleNames(pedigree))
	lrrlist <- lapply(lrrlist, function(x, sampleNames){
		colnames(x) <- sampleNames
		return(x)
	}, sampleNames=sampleNames(pedigree))
	fathers <- paste(fatherNames(pedigree), ext, sep="")
	mothers <- paste(motherNames(pedigree), ext, sep="")
	offsprg <- paste(offspringNames(pedigree), ext, sep="")
	trioindex <- seq_len(nrow(pedigree))
	## want to avoid reading in more than 10 files per node -- would
	## require a lot of total ram, depending on how many nodes are avail.
	if(length(trioindex) > 10){
		ilist <- splitIndicesByLength(trioindex, 10)
	} else{
		ilist <- splitIndicesByNode(trioindex)
		ilist <- ilist[sapply(ilist, function(x) length(x) > 0)]
	}
	if(isPackageLoaded("ff")){
		res <- foreach(i=ilist, .packages="MinimumDistance") %do% MinimumDistance:::read.bsfiles2(path=path, filenames=fathers[i],
						   sampleNames=sampleNames(pedigree)[i],
						   index=index,
						   z=1,
						   baflist=baflist,
						   lrrlist=lrrlist)
		res <- foreach(i=ilist, .packages="MinimumDistance") %do% MinimumDistance:::read.bsfiles2(path=path, filenames=mothers[i],
						   sampleNames=sampleNames(pedigree)[i],
						   index=index,
						   z=2,
						   baflist=baflist,
						   lrrlist=lrrlist)
		res <- foreach(i=ilist, .packages="MinimumDistance") %do% MinimumDistance:::read.bsfiles2(path=path, filenames=offsprg[i],
						   sampleNames=sampleNames(pedigree)[i],
						   index=index,
						   z=3,
						   baflist=baflist,
						   lrrlist=lrrlist)
	} else {
		F <- read.bsfiles2(path=path, filenames=fathers, sampleNames=sampleNames(pedigree))
		M <- read.bsfiles2(path=path, filenames=mothers, sampleNames=sampleNames(pedigree))
		O <- read.bsfiles2(path=path, filenames=offsprg, sampleNames=sampleNames(pedigree))
		for(j in seq_along(index)){
			k <- index[[j]]
			baflist[[j]][, , 1] <- F[k, 2, ]
			baflist[[j]][, , 2] <- M[k, 2, ]
			baflist[[j]][, , 3] <- O[k, 2, ]
			lrrlist[[j]][, , 1] <- F[k, 1, ]
			lrrlist[[j]][, , 2] <- M[k, 1, ]
			lrrlist[[j]][, , 3] <- O[k, 1, ]
		}
	}
	ad <- AssayDataList(logRRatio=lrrlist,
			    BAF=baflist)
	return(ad)
}
