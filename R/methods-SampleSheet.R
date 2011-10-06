SampleSheet <- function(filename,
			id,
			plate, ...){
	df <- data.frame(filename=filename,
			 id=id,
			 plate=plate,
			 ..., stringsAsFactors=FALSE)
	return(as(df, "SampleSheet"))
}
setMethod("sampleNames", signature(object="SampleSheet"),
	  function(object){
		  object$id
	  })
setMethod("show", signature(object="SampleSheet"),
	  function(object){
		  print(head(object))
	  })

setMethod("[", signature(x="SampleSheet"),
	  function(x, i, j, ..., drop=FALSE){
		  x <- callNextMethod()
		  as(x, "SampleSheet")
	  })
