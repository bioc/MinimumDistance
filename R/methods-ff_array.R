## Fix strange behavor for [ with ff_arrays.
setMethod("[", signature(x="ff_array"), function(x, i, j, ..., drop=FALSE){
	if(missing(drop)) drop <- TRUE
	if(length(list(...)) > 0){
		k <- list(...)[[1]]
		if(is(k, "character")) {
			k <- match(dimnames(x)[[3]])
			stopifnot(length(k)>0)
		}
	} else k <- NULL
	if(is.null(k)){
		if(!missing(i) & missing(j)){
			x <- x[i, 1:ncol(x), 1:3, ..., drop=drop]
		}
		if(!missing(i) & !missing(j)){
			x <- x[1:nrow(x), j, 1:3, drop=drop]
		}
		if(missing(i) & !missing(j)){
			x <- x[1:nrow(x), j, 1:3, drop=drop]
		}
	} else {
		if(!missing(i) & missing(j)){
			x <- x[i, 1:ncol(x), k, drop=drop]
		}
		if(!missing(i) & !missing(j)){
			x <- x[i, j, k, drop=drop]
		}
		if(missing(i) & !missing(j)){
			x <- x[1:nrow(x), j, k, drop=drop]
		}
	}
	return(x)
})
