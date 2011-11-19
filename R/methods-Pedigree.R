setMethod("initialize", signature(.Object="Pedigree"),
	  function(.Object, trios=data.frame(),
		   trioIndex=data.frame()){
		  .Object@trios <- trios
		  .Object@trioIndex <- trioIndex
		  return(.Object)
	  })

setValidity("Pedigree", function(object){
	if(!identical(colnames(trios(object)), c("F", "M", "O")))
		return("column names should be 'F', 'M', and 'O'")
	if(any(duplicated(offspringNames(object))))
		return("offspring identifiers must uniquely identify a trio")
	if(any(is.na(unlist(trios(object)))))
		return("Missing values not allowed in pedigree")
	if(!all(c(is(fatherNames(object), "character"),
		  is(motherNames(object), "character"),
		  is(sampleNames(object), "character"))))
		return("sample identifiers must be character strings (e.g., not factors)")
})

setGeneric("trios", function(object) standardGeneric("trios"))
setGeneric("trioIndex", function(object) standardGeneric("trioIndex"))
setMethod("trios", signature(object="Pedigree"),
	  function(object) object@trios)
setMethod("trioIndex", signature(object="Pedigree"),
	  function(object) object@trioIndex)

setMethod("offspringNames", signature(object="Pedigree"), function(object) trios(object)$O)
setMethod("sampleNames", signature(object="Pedigree"), function(object) offspringNames(object))
setMethod("allNames", signature(object="Pedigree"), function(object) trioIndex(object)$individualId)
setMethod("fatherNames", signature(object="Pedigree"), function(object) trios(object)$F)
setMethod("motherNames", signature(object="Pedigree"), function(object) trios(object)$M)
setMethod("show", signature(object="Pedigree"),
	  function(object){
		  cat("trios:\n")
		  print(head(trios(object)))
		  cat("\npedigreeIndex:\n")
		  print(head(trioIndex(object)))
		  cat("\n")
	  })
