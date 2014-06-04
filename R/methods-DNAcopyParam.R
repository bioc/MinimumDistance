DNAcopyParam <- function(alpha=0.01, min.width=2L, undo.splits=c("none", "prune", "sdundo"), undo.SD=3){
  new("DNAcopyParam", alpha=alpha, min.width=min.width, undo.splits=match.arg(undo.splits), undo.SD=undo.SD)
}

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
