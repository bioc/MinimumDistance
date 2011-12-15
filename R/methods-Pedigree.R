setMethod("initialize", signature(.Object="Pedigree"),
	  function(.Object,
		   father=character(),
		   mother=character(),
		   offspring=character(),
		   trios=data.frame(F=father,M=mother,O=offspring),
		   trioIndex=data.frame()){
		  .Object@trios <- trios
		  .Object@trioIndex <- trioIndex
		  return(.Object)
	  })

Pedigree <- function(){
	new("Pedigree")
}

setValidity("Pedigree", function(object){
	if(nrow(trios(object)) > 0){
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
	}
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
			  sns <- unlist(trios(x))
			  x@trioIndex <- trioIndex(x)[match(sns, allNames(x)), ]
		  }
		  return(x)
	  })
