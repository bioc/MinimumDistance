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
