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
