MinDistPosterior <- function(granges=GRangesList()){
  new("MinDistPosterior", granges=granges)
}

setMethod("granges", "MinDistPosterior", function(x, use.mcols=FALSE, ...) x@granges)


setMethod("segs", "MinDistPosterior", function(object){
  unlist(granges(object))
})

setMethod("state", "MinDistPosterior", function(object){
  segs(object)$calls
})

setMethod("cnvFilter", "MinDistPosterior", function(object, filters=FilterParamMD()){
  granges <- segs(object)
  cnvFilter(granges, filters)
})

#' @aliases denovoHemizygous,MinDistPosterior-method
#' @rdname denovo
setMethod("denovoHemizygous", "MinDistPosterior", function(object, filters=FilterParamMD(state="221")){
  cnvFilter(object, filters)
})

#' @param x a \code{MinDistPosterior} object
#' @param i an index for subsetting rows
#' @param j an index for subsetting columns
#' @param ... additional arguments passed to subsetting matrices
#' @param drop logical -- whether to coerce single-row matrices to vectors
#' @aliases [,MinDistPosterior,ANY,ANY,ANY-method
#' @rdname MinDistPosterior-class
setMethod("[", "MinDistPosterior", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@granges <- x@granges[i]
  }
  x
})

##setMethod("[[", "MinDistPosterior", function(x, i, j, ..., drop=FALSE){
##  if(!missing(i)){
##    x@granges <- x@granges[[i]]
##  }
##  x
##})

#' @aliases denovoHomozygous,MinDistPosterior-method
#' @rdname denovo
setMethod("denovoHomozygous", "MinDistPosterior", function(object, filters=FilterParamMD(state="220")){
  cnvFilter(object, filters)
})

#' @aliases denovoDuplication,MinDistPosterior-method
#' @rdname denovo
setMethod("denovoDuplication", "MinDistPosterior", function(object, filters=FilterParamMD(state="224")){
  cnvFilter(object, filters)
})

#' @aliases denovo,MinDistPosterior-method
#' @rdname denovo
setMethod("denovo", "MinDistPosterior", function(object, filters=FilterParamMD(state=c("220", "221", "224"))){
  cnvFilter(object, filters)
})

#' @param object a \code{MinDistPosterior} object
#' @aliases show,MinDistPosterior-method
#' @rdname MinDistPosterior-class
setMethod("show", "MinDistPosterior", function(object){
  cat("Object of class 'MinDistPosterior'\n")
  nsegs <- elementLengths(granges(object))
  cat("  no. segments:\n")
  for(i in seq_along(nsegs)){
    cat("    offspring", i, "-", nsegs[i], "\n")
  }
  g_denovo <- denovo(object)
  denovo_hemi <- denovoHemizygous(object)
  denovo_homo <- denovoHomozygous(object)
  denovo_dup <- denovoDuplication(object)
  cat("  no. denovo     :", length(g_denovo), "\n")
  cat("      hemizygous :", length(denovo_hemi), "\n")
  cat("      homozygous :", length(denovo_homo), "\n")
  cat("      duplication:", length(denovo_dup), "\n")
  cat("See denovoHemizygous(), denovoHomozygous(), denovo()\n")
})
