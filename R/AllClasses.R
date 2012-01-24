setOldClass("ff_array")
setOldClass("ff_matrix")
setClass("LogRratioSet", contains="eSet")
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
#~~ inserted by sgy 1/24/12 ~~~~~~~~~~~~~~~~~~~~~~
setClassUnion("arrayORff_array", c("array", "ff_array"))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setClass("Pedigree", representation(trios="data.frame", trioIndex="data.frame"))
setValidity("Pedigree", function(object){
	msg <- validPedigree(object)
	if(is.null(msg)) return(TRUE) else return(msg)
})

setClass("TrioSet", contains="gSet",
	 representation(fatherPhenoData="AnnotatedDataFrame",
			motherPhenoData="AnnotatedDataFrame",
			pedigree="Pedigree",
			mindist="matrixOrNULL"))

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
		stop("sampleNames of TrioSetList object must be the same as the sampleNames of the phenoData")
	}
	if(!identical(fatherNames(object), sampleNames(fatherPhenoData(object)))){
		stop("fatherNames of TrioSetList object must be the same as the sampleNames of the fatherPhenoData")
	}
	if(!identical(motherNames(object), sampleNames(motherPhenoData(object)))){
		stop("motherNames of TrioSetList object must be the same as the sampleNames of the motherPhenoData")
	}
	if(!is.null(mindist(object))){
		if(!identical(colnames(mindist(object)), sampleNames(object)))
			stop("colnames of mindist matrix must be same as the sampleNames of the TrioSet object")
	}
})

##setMethod("updateObject", signature(object="TrioSetList"),
##          function(object, ..., verbose=FALSE) {
##		  obj <- tryCatch(callNextMethod(), error=function(e) NULL)
##		  if(is.null(obj)){
##			  stop("updateObject failed")
##		  }
##		  return(object)
##	  })


setClass("TrioSetList",
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
	if(length(featureDataList(object)) > 0){
		featureDataClasses <- sapply(featureDataList(object), class)
		if(!unique(featureDataClasses) == "GenomeAnnotatedDataFrame"){
			return("featureDataList must be comprised of GenomeAnnotatedDataFrame(s)")
		}
	}
})



