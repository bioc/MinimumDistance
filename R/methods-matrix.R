setMethod("calculateMindist", signature(object="array"),
	  function(object,..., verbose=TRUE){
		  ##stopifnot(ncol(object)==3)
		  d1 <- object[, , "O"] - object[, , "F"]
		  d2 <- object[, , "O"] - object[, , "M"]
		  I <- as.numeric(abs(d1) <= abs(d2))
		  md <- I*d1 + (1-I)*d2
		  md <- as.matrix(md)
		  return(md)
	  })
