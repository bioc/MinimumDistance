validPedigree <- function(object){
	msg <- NULL
	if(nrow(trios(object)) > 0){
		if(!identical(colnames(trios(object)), c("F", "M", "O"))){
			msg <- "column names should be 'F', 'M', and 'O'"
			return(msg)
		}
		if(any(duplicated(offspringNames(object)))){
			msg <- "offspring identifiers must uniquely identify a trio"
			return(msg)
		}
		if(any(is.na(unlist(trios(object))))){
			msg <- "Missing values not allowed in pedigree"
			return(msg)
		}
		if(!all(c(is(fatherNames(object), "character"),
			  is(motherNames(object), "character"),
			  is(sampleNames(object), "character")))){
			msg <- "sample identifiers must be character strings (e.g., not factors)"
			return(msg)
		}
		if(!all(allNames(object) %in% originalNames(unlist(trios(object))))){
			msg <- "all 'individualId' in slot pedigreeIndex must correspond to an id in the trio slot"
			return(msg)
		}
		if(any(fatherNames(object) == motherNames(object))){
			msg <- "fatherNames can not be the same as the motherNames"
			return(msg)
		}
	}
	return(msg)
}



setMethod("initialize", signature(.Object="Pedigree"),
	  function(.Object, trios, trioIndex, ...){
		  callNextMethod(.Object, trios=trios, trioIndex=trioIndex, ...)
	  })

Pedigree <- function(pedigreeInfo,
		     fatherIds=character(),
		     motherIds=character(),
		     offspringIds=character()){
	if(!missing(pedigreeInfo)){
		msg <- "pedigreeInfo must be a data.frame with column names 'F', 'M', and 'O'"
		if(!is(pedigreeInfo, "data.frame"))
			stop(msg)
		nms <- colnames(pedigreeInfo)
		if(!all(nms %in% c("F", "M", "O"))){
			stop(msg)
		}
		if(any(duplicated(pedigreeInfo$O))){
			stop("offspring identifiers must be unique")
		}
		trios <- data.frame(F=make.unique2(as.character(pedigreeInfo$F)),
				    M=make.unique2(as.character(pedigreeInfo$M)),
				    O=as.character(pedigreeInfo$O),
				    stringsAsFactors=FALSE)
		allIds <- as.character(unlist(trios))
	} else {
		if(!(length(fatherIds)==length(motherIds) & length(fatherIds)==length(offspringIds))){
			stop("father, mother, and offspring identifiers must be the same length")
		}
		fatherIds <- as.character(fatherIds)
		motherIds <- as.character(motherIds)
		offspringIds <- as.character(offspringIds)
		trios <- data.frame(F=make.unique2(fatherIds),
				    M=make.unique2(motherIds),
				    O=offspringIds,
				    stringsAsFactors=FALSE)
		allIds <- c(fatherIds, motherIds, offspringIds)
		if(any(duplicated(offspringIds))){
			stop("Can not have duplicated offspring identifiers. An offspring can only belong to one trio.")
		}
	}
	trio.index <- as.integer(matrix(seq_len(nrow(trios)), nrow(trios), 3, byrow=FALSE))
	memberId <- rep(c("F", "M", "O"), each=nrow(trios))
	pedigreeIndex <- data.frame(individualId=allIds,
				    memberId=memberId,
				    index.in.pedigree=trio.index,
			    stringsAsFactors=FALSE)
	rownames(pedigreeIndex) <- NULL
	new("Pedigree", trios=trios, trioIndex=pedigreeIndex)
}


setMethod("trios", signature(object="Pedigree"),
	  function(object) object@trios)
setMethod("trioIndex", signature(object="Pedigree"),
	  function(object) object@trioIndex)
setMethod("offspringNames", signature(object="Pedigree"), function(object) trios(object)$O)
setMethod("sampleNames", signature(object="Pedigree"), function(object) offspringNames(object))
setMethod("allNames", signature(object="Pedigree"), function(object) unique(trioIndex(object)$individualId))
setMethod("fatherNames", signature(object="Pedigree"), function(object) trios(object)$F)
setMethod("motherNames", signature(object="Pedigree"), function(object) trios(object)$M)
setMethod("show", signature(object="Pedigree"),
	  function(object){
		  cat("trios:\n")
		  print(head(trios(object)))
		  if(nrow(trios(object)) > 6)
			  cat(".\n.\n.\n")
		  cat("\npedigreeIndex:\n")
		  print(head(trioIndex(object)))
		  if(nrow(trioIndex(object)) > 6)
			  cat(".\n.\n.")
		  cat("\n")
	  })

setMethod("[", signature(x="Pedigree"),
	  function(x, i, j, ..., drop=FALSE){
		  if(missing(i)) {
			  return(x)
		  } else {
			  x@trios <- trios(x)[i, ]
			  sns <- originalNames(unlist(trios(x)))
			  x@trioIndex <- trioIndex(x)[match(sns, trioIndex(x)$individualId), ]
		  }
		  return(x)
	  })

setMethod("dim", signature(x="Pedigree"), function(x){
	dim(trios(x))
})

setMethod("annotatedDataFrameFrom", signature(object="Pedigree", byrow="logical"),
	  function (object, byrow, sample.sheet, which=c("offspring", "father", "mother"),
		    row.names=NULL, ...){
		  dims <- dim(object)
		  if (is.null(dims) || all(dims == 0)){
			  return(annotatedDataFrameFrom(NULL, byrow = byrow, ...))
		  }
		  which <- match.arg(which)
		  nms <- switch(which,
				offspring=offspringNames(object),
				father=fatherNames(object),
				mother=motherNames(object))
		  nms <- make.unique2(nms)
		  if(missing(sample.sheet)){
			  n <- length(nms)
			  data <- data.frame(numeric(n), row.names = nms)[, FALSE]
			  dimLabels <-  c("sampleNames", "sampleColumns")
			  phenoData <- new("AnnotatedDataFrame", data = data, dimLabels = dimLabels)
		  } else {
			  if(is.null(row.names)){
				  stop("sample.sheet is not missing, but row.names not specified")
			  } else{
				  stopifnot(originalNames(nms) %in% row.names)
				  index <- match(originalNames(nms), row.names)
				  data <- sample.sheet[index, ]
				  rownames(data) <- nms
				  dimLabels <-  c("sampleNames", "sampleColumns")
				  phenoData <- new("AnnotatedDataFrame", data = data, dimLabels = dimLabels)
			  }
		  }
		  return(phenoData)
	  })




