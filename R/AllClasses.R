setOldClass("ff_array")
setOldClass("ff_matrix")
setClassUnion("matrixOrNULL", c("matrix", "NULL", "ff_matrix"))
setClassUnion("matrixOrff", c("matrix", "ff_matrix"))
setClassUnion("arrayOrNULL", c("array", "NULL"))
setClass("LogRratioSet", contains="eSet")
setClass("LikSet",
	 contains="LogRratioSet",
	 representation(loglik="array",
			range.index="integer"),
	 prototype = prototype(
	 new("VersionedBiobase",
	     versions=c(classVersion("LogRratioSet"), LikSet="1.0.0"))))


##setClass("SampleSheet", contains="DataFrame")


setClass("Pedigree", representation(trios="data.frame", trioIndex="data.frame"))

setValidity("Pedigree", function(object){
	msg <- validPedigree(object)
	if(is.null(msg)) return(TRUE) else return(msg)
})

##setClass("SampleSheet", representation(phenoData="list"))
##setValidity("SampleSheet", function(object){
##	if(length(phenoData(object)) != 3){
##		return("length of SampleSheet object should be 3 (one DataFrame for each member in the trio)")
##	}
##	if(!all(sapply(phenoData(object), class) == "DataFrame")){
##		return("Each element of the SampleSheet must be a DataFrame")
##	}
##	d <- lapply(phenoData(object), dim)
##	firstElement <- d[[1]]
##	d <- d[-1]
##	if(!all(sapply(d, function(x) x == firstElement))){
##		return("Each element of the SampleSheet list must have the same dimension")
##	}
##	nms <- lapply(phenoData(object), rownames)
##	if(!all(sapply(nms, function(x) identical(x, nms[[1]])))){
##		return("rownames must be the same for each element in phenoData")
##	}
##})

##setClass("TrioSet", contains="LogRratioSet",
##	 representation(phenoData2="arrayOrNULL",
##			mindist="matrixOrNULL",
##			mad="matrix"),
##	 prototype = prototype(
##	                       new("VersionedBiobase",
##				   versions=c(classVersion("LogRratioSet"), TrioSet="0.0.4"))))
##
##setClass("TrioSet", contains="LogRratioSet", ##contains="LogRratioSet",
##	 representation(phenoData2="array",
##			mindist="matrixOrNULL",
##			mad="matrix"),
##	 prototype = prototype(
##	                       new("VersionedBiobase",
##				   versions=c(classVersion("LogRratioSet"), TrioSet="0.0.5"))))

##setClass("TrioSet", contains="LogRratioSet", ##contains="LogRratioSet",
setClass("TrioSet", contains="eSet",
	 representation(##sampleSheet="SampleSheet",
			fatherPhenoData="AnnotatedDataFrame",
			motherPhenoData="AnnotatedDataFrame",
			pedigree="Pedigree"))

setValidity("TrioSet", function(object){
	ped <- pedigree(object)
	##ss <- sampleSheet(object)
	validObject(ped)
	##validObject(ss)
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
		stop("sampleNames of TrioSetList object must be the same as the sampleNames of the phenoData")
	}
	if(!identical(fatherNames(object), sampleNames(fatherPhenoData(object)))){
		stop("fatherNames of TrioSetList object must be the same as the sampleNames of the fatherPhenoData")
	}
	if(!identical(motherNames(object), sampleNames(motherPhenoData(object)))){
		stop("motherNames of TrioSetList object must be the same as the sampleNames of the motherPhenoData")
	}
})

##setMethod("updateObject", signature(object="TrioSet"),
##          function(object, ..., verbose=FALSE) {
##		  obj <- tryCatch(callNextMethod(), error=function(e) NULL)
##		  if(is.null(obj)){
##			  stop("updateObject failed")
##		  }
##		  return(object)
##	  })

## TrioSetList has assay data elements that are lists of arrays
## TrioSet has assay data elements that are arrays.
##
## object[[1]] returns a TrioSet
## object[c(1, 2)] returns a TrioSetList
## object[c(1, 1)] returns a TrioSetList
## The featureData is also a list.
## Require that all elements of the list have the same number of columns
##
## Redefine TrioSet class such that an element of the assayData is either a list or an array.
##
## Each element in the list is the assay data for a single chromosome.
## the accessor lrr() returns a list
## the lrr()[[chr]] gets the info for one chromosome
##
## object[[1]]
##

##setClass("TrioSetList", contains="list",
##	 representation(pedigree="Pedigree",
##			sampleSheet="SampleSheet"))
## Do away with SampleSheet. Just use phenoData...
##setClass("AssayDataList", contains="AssayData")
##setClassUnion("AssayDataList", c("list", "environment"))
##setMethod("dim", "AssayDataList", function(x) c(length(x), 0))


##AssayDataList <- function(storage.mode = c("lockedEnvironment", "environment", "list"), ...) {
##	storage.mode <- match.arg(storage.mode) ## defaults to "lockedEnvironment"
##	assayData <- switch(storage.mode,
##			    lockedEnvironment =,
##			    environment = new.env(parent=emptyenv()),
##			    list = list())
##	arglist <- list(...)
##	for (nm in names(arglist)) assayData[[nm]] <- arglist[[nm]]
##	if (storage.mode == "lockedEnvironment") Biobase:::assayDataEnvLock(assayData)
##	assayData
##}

setClass("TrioSetList", ##contains="eSet",
	 representation(pedigree="Pedigree",
			assayDataList="AssayData",
			phenoData="AnnotatedDataFrame",
			fatherPhenoData="AnnotatedDataFrame",
			motherPhenoData="AnnotatedDataFrame",
			featureDataList="list",
			chromosome="integer"))

setValidity("TrioSetList", function(object){
	nms <- ls(assayData(object))
	if(!all(c("BAF", "logRRatio") %in% nms)){
		msg <- "BAF and logRRatio are required elements of the assayData"
		return(msg)
	}
	if(length(object) > 0){
		msg <- validAssayDataDims(assayData(object))
		if(!all(msg == TRUE)) return(msg)
		elt <- (ls(assayDataList(object)))[[1]]
		b <- assayDataList(object)[[elt]]
		if(length(chromosome(object)) != length(b)){
			return("chromosome slot must be the same length as the length of the list for each assayData element")
		}
	}
	validObject(pedigree(object))
	if(!identical(sampleNames(object), sampleNames(phenoData(object)))){
		stop("sampleNames of TrioSetList object must be the same as the sampleNames of the phenoData")
	}
	if(!identical(fatherNames(object), sampleNames(fatherPhenoData(object)))){
		stop("fatherNames of TrioSetList object must be the same as the sampleNames of the fatherPhenoData")
	}
	if(!identical(motherNames(object), sampleNames(motherPhenoData(object)))){
		stop("motherNames of TrioSetList object must be the same as the sampleNames of the motherPhenoData")
	}
	if(length(featureDataList(object)) != length(chromosome(object))){
		return("each chromosome should have an element in the featureDataList")
	}
})


##setMethod("dim", "AssayDataList", function(x) c(length(x), 0))

##setClass("TrioSetList", contains="eSet", representation(assayData="AssayDataList"))


##tmp=new("SimpleList")
##trace(assayDataValidMembers, browser)
##ad <- assayDataNew(logRRatio=SimpleList(), BAF=SimpleList())
##ad <- assayDataNew(logRRatio=matrix(), BAF=matrix())
##tmp=new("TrioSetList")

##setValidity("TrioSetList", function(object){
##	if(!all(unlist(pedigree(object)) %in% sampleNames(sampleSheet(object))))
##		return("All names in the pedigree object must be present in the sample sheet")
##	if(length(object) > 0){
##		if(!identical(sampleNames(pedigree(object)), colnames(lrr(object[[1]]))))
##			return("The sampleNames of the pedigree slot must be the same as the column names of the assayData elements in the TrioSet")
##	}
##	if(!identical(allNames(pedigree(object)), sampleNames(sampleSheet(object))))
##		return("allNames of Pedigree must be identical to sampleNames of SampleSheet")
##	if(!identical(as.character(unlist(trios(object))), as.character(phenoData2(object[[1]])[, "sampleNames", ])))
##		return("The phenoData2 slot for the elements in the TrioSetList must have a sampleNames column equal to trios(object)")
##	if(!checkOrder(object)) return("each element in the TrioSetList must be ordered by chromosome and physical position.")
##	TRUE
##})


