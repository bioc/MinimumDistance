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
