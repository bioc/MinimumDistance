setMethod("calculateMindist", signature(object="list"),
	  function(object){
		  lapply(object, calculateMindist)
	  })
