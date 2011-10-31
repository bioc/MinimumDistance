##RangedDataMinimumDistance <- function(ranges=IRanges(), ...){
##	browser()
##	rd <- RangedDataCNV(ranges=ranges, ...)
##	new("RangedDataMinimumDistance", ranges=ranges(rd), values=values(rd))
##}
##RangedDataCBS2 <- function(ranges=IRanges(),
##			   state=vector("character", length(ranges)), ...){
##	rd <- RangedDataCBS(ranges=ranges, state=state, ...)
##	new("RangedDataCBS2", ranges=ranges(rd), values=values(rd))
##}
##setMethod("trioNames", signature(object="RangedDataCNV"), function(object) {
##	res <- object$family
##	if(is.null(res)) stop("'family' not in the column names")
##	return(res)
##})
##setMethod("RangedDataCNV", signature(ranges="IRanges"),
##	  function(ranges=IRanges(), ...,
##		   space=NULL,
##		   universe=NULL){
##		  nms <- names(list(...))
##		  stopifnot(c("chrom", "id", "num.mark", "seg.mean", "start.index", "end.index") %in% nms)
##		  rd <- RangedData(ranges=ranges, ..., space=space, universe=universe)
##		  rd2 <- as(rd, "RangedDataCNV")
##		  return(rd2)
##	  })

##setMethod("plot", signature(x="RangedDataCNV", y="missing"),
##	  function(x, y, ...){
##		  df <- todf(x)
##		  plot(x=df, ...)
##	  })


##setMethod("[", signature(x="RangedDataCNV"),
##	  function(x, i, j, ..., drop=FALSE){
##		  ## The "[" method for RangedData does not have separate methods for the slot elements
##		  ## It seems that the easiest way to subset an object extending the RangedData class is to use
##		  ## the RangeData method directly -- this will return an object of class RangedData
##		  rd <- callNextMethod(x, i, j, ..., drop=drop)
##		  ## Now, update the components of 'x'
##		  x@ranges <- ranges(rd)
##		  x@values <- values(rd)
##		  return(x)
##	  })



