setOldClass("ff_array")
setOldClass("ff_matrix")
setClassUnion("matrixOrNULL", c("matrix", "NULL", "ff_matrix"))
setClassUnion("matrixOrff", c("matrix", "ff_matrix"))
setClassUnion("arrayOrNULL", c("array", "NULL"))
setClass("LogRratioSet", contains="eSet")
##setValidity("LogRratioSet", function(object){
##	#assayDataValidMembers(assayData(object), c("logRRatio", "BAF"))
##})
setClass("LikSet",
	 contains="LogRratioSet",
	 representation(loglik="array",
			range.index="integer"),
	 prototype = prototype(
	 new("VersionedBiobase",
	     versions=c(classVersion("LogRratioSet"), LikSet="1.0.0"))))
setClass("TrioSet", contains="LogRratioSet",
	 representation(phenoData2="arrayOrNULL",
			mindist="matrixOrNULL",
			mad="matrix"),
	 prototype = prototype(
	                       new("VersionedBiobase",
				   versions=c(classVersion("LogRratioSet"), TrioSet="0.0.4"))))

setClass("TrioSet", contains="LogRratioSet", ##contains="LogRratioSet",
	 representation(phenoData2="array",
			mindist="matrixOrNULL",
			mad="matrix"),
	 prototype = prototype(
	                       new("VersionedBiobase",
				   versions=c(classVersion("LogRratioSet"), TrioSet="0.0.5"))))
## should we add a slot for trioNames
##  -- would be a R x 3 matrix, where R is the number of trios
##  -- R must be equal to the number of columns of the assayData arrays
##  -- '[' method for the trioNames slot
##  -- 'show' method for trioNames slot
##  -- add trioNames accessor / replacement method
##  -- modify offspringNames, fmoNames, ... to access trioNames slot
## Other potential slots:
##  -- dna source
##  -- batch / plate

setMethod("updateObject", signature(object="TrioSet"),
          function(object, ..., verbose=FALSE) {
		  obj <- tryCatch(callNextMethod(), error=function(e) NULL)
		  if(is.null(obj)){
			  stop("updateObject failed")
		  }
		  return(object)
	  })
setClass("SampleSheet", contains="DataFrame")

## might try extending AnnotatedDataFrame instead of data.frame
setClass("Pedigree", contains="list",
	 representation(trios="data.frame",
			trioIndex="data.frame"))





setClass("TrioSetList", contains="list",
	 representation(pedigree="Pedigree",
			sampleSheet="SampleSheet"))

setValidity("TrioSetList", function(object){
	if(!all(unlist(pedigree(object)) %in% sampleNames(sampleSheet(object))))
		return("All names in the pedigree object must be present in the sample sheet")
	if(length(object) > 0){
		if(!identical(sampleNames(pedigree(object)), colnames(lrr(object[[1]]))))
			return("The sampleNames of the pedigree slot must be the same as the column names of the assayData elements in the TrioSet")
	}
	if(!identical(allNames(pedigree(object)), sampleNames(sampleSheet(object))))
		return("allNames of Pedigree must be identical to sampleNames of SampleSheet")
	if(!identical(as.character(unlist(trios(object))), as.character(phenoData2(object[[1]])[, "sampleNames", ])))
		return("The phenoData2 slot for the elements in the TrioSetList must have a sampleNames column equal to trios(object)")
	if(!checkOrder(object)) return("each element in the TrioSetList must be ordered by chromosome and physical position.")
	TRUE
})


