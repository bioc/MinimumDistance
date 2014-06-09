#' @export
MinDistParam <- function(nMAD=0.75, dnacopy=DNAcopyParam(), penncnv=PennParam()){
  new("MinDistParam", nMAD=nMAD, dnacopy=dnacopy, penncnv=penncnv)
}


setMethod("nMAD", "MinDistParam", function(object) object@nMAD)

setReplaceMethod("nMAD", c("MinDistParam", "numeric"), function(object, value){
  object@nMAD <- value
  object
})

setMethod("penncnv", "MinDistParam", function(object) object@penncnv)

##setGeneric("penncnv", function(object) standardGeneric("penncnv"))
##setMethod("penncnv", "MinDistParam", function(object) object@penncnv)
##setReplaceMethod("penncnv", c("MinDistParam", "PennParam"),
##                 function(object, value){
##                   object@penncnv <- value
##                   object
##                 })

setMethod("show", "MinDistParam", function(object){
  cat("An object of class 'MinDistParam'\n")
  cat("  call segments with |seg.mean|/MAD > nMAD\n")
  cat("  nMAD = ", nMAD(object), "\n")
  cat("DNAcopy params:\n")
  p <- dnacopy(object)
  cat("    alpha: ", alpha(p), "\n")
  cat("    min.width: ", min.width(p), "\n")
  cat("    undo.splits: ", undo.splits(p), "\n")
  cat("    undo.SD: ", undo.SD(p), "\n")
  cat(" Setting nMAD() to smaller values will increase the number of segments that are called.\n")
  cat(" See ?segment for description of DNAcopy parameters\n")
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
