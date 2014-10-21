#' Parameters for filtering results from the segmentation and copy number inference
#'
#' A container for criteria used to filter the segmentation results
#' post-hoc.  Options including filtering on the posterior call, the
#' posterior probability of the posterior call, the minimum number of
#' markers spanned by the segment, the minimum width of the segment,
#' and chromosome.  Convenience functions are available for commonly
#' used filters.
#' @param state trio copy number states to select
#' @param seqnames chromosome names to select
#' @param ... additional arguments passed to \code{\link[VanillaICE]{FilterParam}}
#' @examples
#' library(VanillaICE)
#' data(md_gr)
#' data(md_exp)
#' mdparam <- MinDistParam()
#' fit <- MAP2(md_exp, md_gr, mdparam)
#' ## return all segments
#' segs(fit)
#'
#' ## Default filters
#' param <- FilterParamMD()
#' param
#' cnvFilter(fit, param)
#'
#' param2 <- FilterParamMD(seqnames="chr22", probability=0.9, numberFeatures=10)
#' cnvFilter(fit, param2)
#' denovoHemizygous(fit)
#' @rdname FilterParamMD
#' @export
FilterParamMD <- function(state=trioStateNames(), seqnames=paste0("chr", 1:22), ...){
  param <- FilterParam(state=state, seqnames=seqnames, ...)
  as(param, "FilterParamMD")
}


#' @param object a \code{FilterParamMD} object
#' @aliases show,FilterParamMD-method
#' @rdname FilterParamMD-class
#' @export
setMethod("show", "FilterParamMD", function(object){
  cat("An object of class 'FilterParamMD'\n")
  cat("   min. posterior probability of trio CNV call:", probability(object), "\n")
  cat("   min. no. of markers spanned by segment     :", numberFeatures(object), "\n")
  cat("   min. width of segment                      :", width(object), "\n")
  states <- state(object)
  nstates <- length(states)
  if(length(states) > 6){
    states <- paste0(c(states[1:6], "..."), collapse=", ")
  } else states <- paste0(states, collapse=", ")
  cat("  ",nstates, "selected trio CN states              :", states, "\n")
  chroms <- chromosome(object)
  if(length(chroms) > 6){
    chroms <- paste(paste(chroms[1:6], collapse=", "), "...")
  } else paste(chroms, collapse=", ")
  cat("   selected seqnames                         :", chroms, "\n")
})
