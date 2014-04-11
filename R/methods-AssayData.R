AssayDataList <- function(storage.mode = c("lockedEnvironment", "environment", "list"), ...) {
  storage.mode <- match.arg(storage.mode) ## defaults to "lockedEnvironment"
  assayData <- switch(storage.mode,
                      lockedEnvironment =,
                      environment = new.env(parent=emptyenv()),
                      list = list())
  arglist <- list(...)
  for (nm in names(arglist)) assayData[[nm]] <- arglist[[nm]]
  ##if (storage.mode == "lockedEnvironment") Biobase:::assayDataEnvLock(assayData)
  storageMode(assayData) <- storage.mode
  assayData
}

assayDataListDims <- function(object) {
  nms <- if (storageMode(object) == "list") names(object) else ls(object)
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

assayDataListLD <- function(path="", ext="", pedigree, featureData, ffprefix=""){
	filenames <- paste(originalNames(allNames(pedigree)), ext, sep="")
	fnames <- file.path(path, filenames)
	stopifnot(all(file.exists(fnames)))
	if(missing(featureData)) stop("featureData can not be missing")
	if(missing(pedigree)) stop("pedigree can not be missing")
	index <- split(seq_len(nrow(featureData)), chromosome(featureData))
	useff <- isPackageLoaded("ff")
	if(useff){
		message("Initializing ff arrays for BAFs and LRRs  ...")
	} else {
		message("ff package not loaded -- initializing arrays for BAFs and LRRs ...")
	}
	pkgs <- if(useff) c("ff", neededPkgs()) else neededPkgs()
	outdir <- ldPath()
	i <- NULL
	bafAndLrrList <- foreach(i=index, .packages=pkgs) %dopar% {
          initializeLrrAndBafArrays(dims=c(length(i), nrow(pedigree), 3), outdir=outdir,
                                    col.names=sampleNames(pedigree),
                                    name=ffprefix)
	}
	baflist <- lapply(bafAndLrrList, "[[", 1)
	lrrlist <- lapply(bafAndLrrList, "[[", 2)
##	} else {
##		baflist <- lapply(index, function(x) initializeBigArray("baf", dim=c(length(x), nrow(pedigree), 3), vmode="integer"))
##		lrrlist <- lapply(index, function(x) initializeBigArray("lrr", dim=c(length(x), nrow(pedigree), 3), vmode="integer"))
##		baflist <- lapply(baflist, function(x, sampleNames) {
##			colnames(x) <- sampleNames
##			return(x)
##		}, sampleNames=sampleNames(pedigree))
##		lrrlist <- lapply(lrrlist, function(x, sampleNames){
##			colnames(x) <- sampleNames
##			return(x)
##		}, sampleNames=sampleNames(pedigree))
##	}
	if(useff)
		message("Reading ", length(fnames),
			" files and writing ff files to directory ", ldPath())
	fathers <- paste(fatherNames(pedigree), ext, sep="")
	mothers <- paste(motherNames(pedigree), ext, sep="")
	offsprg <- paste(offspringNames(pedigree), ext, sep="")
	trioindex <- seq_len(nrow(pedigree))
	## want to avoid reading in more than 10 files per node -- would
	## require a lot of total ram, depending on how many nodes are avail.
	##
	## for reading in data, we can't split by chromosome (all
	## markers are read in at once) So, we split by samples.
	if(parStatus()){
		## e.g., for 900 trios and 3 workers,
		## each worker reads in 300 trios
		##  --- below we read 100 fathers, 100 mothers, 100 offspring...
		##ilist <- splitIndicesByLength(trioindex, getCluster())
		ilist <- splitIndicesByNode(trioindex)
		ilist <- ilist[sapply(ilist, function(x) length(x) > 0)]
	} else {
		## execution is sequential.
		ilist <- list(trioindex)
	}
	if(useff){
		## pass the ff object / array to each worker
		## read in the files and assign the results to column z
		## workers read in different sets of files and assign to the baflist and lrrlist ff objects
		## read one file.
		if(!getDoParWorkers()) registerDoSEQ()
		fns <- featureNames(featureData)
		res <- foreach(i=ilist, .packages=pkgs) %dopar% {
                  read.bsfiles2(path=path,
                                filenames=originalNames(fathers[i]),
                                sampleNames=sampleNames(pedigree)[i],
                                marker.index=index,
                                z=1,
                                baflist=baflist,
                                lrrlist=lrrlist,
                                featureNames=fns)
		}
		res <- foreach(i=ilist, .packages=pkgs) %dopar% {
                  read.bsfiles2(path=path,
                                filenames=originalNames(mothers[i]),
                                sampleNames=sampleNames(pedigree)[i],
                                marker.index=index,
                                z=2,
                                baflist=baflist,
                                lrrlist=lrrlist,
                                featureNames=fns)
		}
		res <- foreach(i=ilist, .packages=pkgs) %dopar% {
			read.bsfiles2(path=path,
                                      filenames=offsprg[i],
                                      sampleNames=sampleNames(pedigree)[i],
                                      marker.index=index,
                                      z=3,
                                      baflist=baflist,
                                      lrrlist=lrrlist,
                                      featureNames=fns)
                      }
		message("Finished reading/writing processed data.")
		gc()
	} else {
		F <- read.bsfiles2(path=path, filenames=originalNames(fathers), sampleNames=sampleNames(pedigree))
		M <- read.bsfiles2(path=path, filenames=originalNames(mothers), sampleNames=sampleNames(pedigree))
		O <- read.bsfiles2(path=path, filenames=offsprg, sampleNames=sampleNames(pedigree))
		## the featureData is ordered by chromosome and physical position
##		.i <- rownames(F) %in% featureNames(featureData)
##		F <- F[.i,,,drop=FALSE]
##		M <- M[.i,,,drop=FALSE]
##		O <- O[.i,,,drop=FALSE]
		mindex <- match(featureNames(featureData), rownames(F))
		F <- F[mindex, , , drop=FALSE]
		M <- M[mindex, , , drop=FALSE]
		O <- O[mindex, , , drop=FALSE]
		##featureData <- featureData[mindex, ]
		##if(!identical(featureNames(featureData), rownames(F))) stop("feature names in featureData do not match rownames in assay data array")
		##index <- split(seq_len(nrow(featureData)), chromosome(featureData))
		for(j in seq_along(index)){
			k <- index[[j]]
			baflist[[j]][, , 1] <- integerMatrix(F[k, 2, ], scale=1000)
			baflist[[j]][, , 2] <- integerMatrix(M[k, 2, ], scale=1000)
			baflist[[j]][, , 3] <- integerMatrix(O[k, 2, ], scale=1000)
			lrrlist[[j]][, , 1] <- integerMatrix(F[k, 1, ], scale=100)
			lrrlist[[j]][, , 2] <- integerMatrix(M[k, 1, ], scale=100)
			lrrlist[[j]][, , 3] <- integerMatrix(O[k, 1, ], scale=100)
		}
	}
	ad <- AssayDataList(logRRatio=lrrlist,
			    BAF=baflist)
	return(ad)
}
