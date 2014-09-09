#' Constructor for \code{MinDistParam} class
#'
#' The \code{MinDistParam} class contains parameters used for the
#' segmentation implemented in the \code{DNAcopy} package, parameters
#' extracted from the \code{PennCNV} HMM such as parent-offspring
#' transmission parobabilities (see citation below), and initial
#' values / parameters for computing emission probabilities.
#'
#' @param nMAD a length-one numeric vector indicating the minimal
#' number of median absolute deviations of the mean segmented minimum
#' distance from zero.  For non-zero segments (# median absolute
#' deviations > nMAD), maximum a posteriori estimates of the
#' parent-offspring copy number states are computed.  Segments with
#' minimum distance values near zero are not called as they are less
#' likely to correspond to regions with de novo copy number
#' alterations.
#' @param dnacopy an object of class \code{DNAcopyParam}.
#' @param penncnv probabilities/parameters of the PennCNV hidden
#' Markov model
#' @param emission an object of class \code{EmissionParam}
#' @param thin a length-one vector indicating whether to thin the
#' data. This is primarily for internal use in conjunction with the
#' \code{filterExperiment} function.
#' @export
MinDistParam <- function(nMAD=0.75, dnacopy=DNAcopyParam(), penncnv=PennParam(), emission=EmissionParam(), thin=10L){
  new("MinDistParam", nMAD=nMAD, dnacopy=dnacopy, penncnv=penncnv, emission=emission, thin=thin)
}

setMethod("thin", "MinDistParam", function(object) object@thin)
setMethod("nMAD", "MinDistParam", function(object) object@nMAD)

setReplaceMethod("nMAD", c("MinDistParam", "numeric"), function(object, value){
  object@nMAD <- value
  object
})

setMethod("stateNames", "MinDistParam", function(object) stateNames(penncnv(object)))

setMethod("penncnv", "MinDistParam", function(object) object@penncnv)
setMethod("emission", "MinDistParam", function(object) object@emission)
setReplaceMethod("emission", c("MinDistParam", "EmissionParam"),
                 function(object, value) {
                   object@emission <- value
                   object
                 })
setMethod("EMupdates", "MinDistParam", function(object) EMupdates(object))

##setGeneric("penncnv", function(object) standardGeneric("penncnv"))
##setMethod("penncnv", "MinDistParam", function(object) object@penncnv)
##setReplaceMethod("penncnv", c("MinDistParam", "PennParam"),
##                 function(object, value){
##                   object@penncnv <- value
##                   object
##                 })

setMethod("show", "MinDistParam", function(object){
  cat("An object of class 'MinDistParam'\n")
  cat("  call segments with |seg.mean|/MAD > ", nMAD(object), "\n")
  cat("  Setting nMAD() to smaller values will increase the number of segments that are called.\n")
  cat("  DNAcopy settings:\n")
  p <- dnacopy(object)
  cat("    alpha: ", alpha(p), "\n")
  cat("    min.width: ", min.width(p), "\n")
  cat("    undo.splits: ", undo.splits(p), "\n")
  cat("    undo.SD: ", undo.SD(p), "\n")
  cat("    See segment() for description of DNAcopy parameters\n")
})

DNAcopyParam <- function(alpha=0.01, min.width=2L, undo.splits=c("none", "prune", "sdundo"), undo.SD=3){
  new("DNAcopyParam", alpha=alpha, min.width=min.width, undo.splits=match.arg(undo.splits), undo.SD=undo.SD)
}

dnacopy <- function(object) object@dnacopy

setMethod("alpha", "DNAcopyParam", function(object) object@alpha)
setMethod("min.width", "DNAcopyParam", function(object) object@min.width)
setMethod("undo.splits", "DNAcopyParam", function(object) object@undo.splits)
setMethod("undo.SD", "DNAcopyParam", function(object) object@undo.SD)

setMethod("show", "DNAcopyParam", function(object){
  cat("An object of class 'DNAcopyParam'\n")
  cat("  alpha: ", alpha(object), "\n")
  cat("  min.width: ", min.width(object), "\n")
  cat("  undo.splits: ", undo.splits(object), "\n")
  cat("  undo.SD: ", undo.SD(object), "\n")
})
