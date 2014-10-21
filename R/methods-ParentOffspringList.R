setMethod("mdgrFile", "ParentOffspringList", function(object, label="mdgr"){
  fileName(object, label)
})

setMethod("mapFile", "ParentOffspringList", function(object, label="map"){
  fileName(object, label)
})

setMethod("sapply", "ParentOffspringList", function(X, FUN, ...,
                                                    simplify=TRUE, USE.NAMES=TRUE){
  FUN <- match.fun(FUN)
  answer <- lapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer)))
    names(answer) <- X
  if (!identical(simplify, FALSE) && length(answer))
    simplify2array(answer, higher = (simplify == "array"))
  else answer
})

setMethod("lapply", "ParentOffspringList", function(X, FUN, ...){
  FUN <- match.fun(FUN)
  J <- seq_len(length(X))
  j <- NULL
  answer <- foreach(j = J, .packages=c("VanillaICE")) %dopar% {
    FUN(X[[j]], ...)
  }
  answer
})

#' Constructor for ParentOffspringList class
#'
#' @param pedigrees a list of \code{ParentOffspring} objects
#' @param id identifier for a pedigree
#' @rdname ParentOffspringList-class
#' @export
ParentOffspringList <- function(pedigrees=list(), id){
  if(missing(id))
    id <- as.character(sapply(pedigrees, pedigreeName))
  new("ParentOffspringList", id=id, pedigrees=pedigrees)
}

#' @param object a \code{ParentOffspringList} object
#' @aliases pedigreeName,ParentOffspringList-method
#' @rdname ParentOffspringList-class
setMethod("pedigreeName", "ParentOffspringList", function(object) object@id)

#' @aliases show,ParentOffspringList-method
#' @rdname ParentOffspringList-class
setMethod("show", "ParentOffspringList", function(object){
  ids <- paste(head(pedigreeName(object)), collapse=", ")
  if(length(pedigreeName(object)) > 6) ids <- paste0(ids, ",...")
  cat("# pedigrees:", length(pedigreeName(object)), "\n")
})

#' @param x a \code{ParentOffspringList} object
#' @param i a numeric vector for subsetting the list (optional)
#' @param j ignored
#' @param ... ignored
#' @param drop ignored
#' @aliases "[[",ParentOffspringList,ANY,ANY-method
#' @rdname ParentOffspringList-class
setMethod("[[", "ParentOffspringList", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x <- x@pedigrees[[i]]
  }
  x
})

#' @aliases [,ParentOffspringList,ANY-method
#' @rdname ParentOffspringList-class
setMethod("[", "ParentOffspringList", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@pedigrees <- x@pedigrees[i]
    x@id <- x@id[i]
  }
  x
})

#' @aliases length,ParentOffspringList-method
#' @rdname ParentOffspringList-class
setMethod("length", "ParentOffspringList", function(x) length(x@id))


setMethod("fileName", "ParentOffspringList", function(object, label){
  datadir <- sapply(object@pedigrees, function(x) dirname(parsedPath(x)[1]), USE.NAMES=FALSE)
  file.path(datadir, paste0(pedigreeName(object), "_", label, ".rds"))
})
