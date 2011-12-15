#setValidity("SampleSheet", function(object){
##	if(nrow(object) > 0){
##		if(!"id" %in% colnames(object))
##			return("'id' needs to be in the column name")
##		if(any(duplicated(object$id)))
##			return("'id' needs to be unique")
##		NULL
##	}
#})
SampleSheet <- function(...){
	dF <- DataFrame(...)
	as(dF, "SampleSheet")
}
setMethod("sampleNames", signature(object="SampleSheet"),
	  function(object){
		  object@rownames
	  })
setMethod("show", signature(object="SampleSheet"),
	  function(object){
		  cat(head(sampleNames(object)))
		  if(length(sampleNames(object)) > 0){
			  cat("...\n")
		  }
	  })

