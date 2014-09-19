
##setGeneric("ParentOffspring", function(object, id=character(),
##                                       father=character(),
##                                       mother=character(),
##                                       offspring=character(),
##                                       parsedPath=character())
##           setGeneric("ParentOffspring"))


#' Constructor for ParentOffspring class
#'
#' @param id length-one character vector providing a family-level id
#' @param father length-one character vector providing sample ids for father
#' @param mother length-one character vector providing sample ids for mother
#' @param offspring  character vector providing sample ids for offspring (can have length greater than one if there is more than one offspring)
#' @param parsedPath character vector providing path to low-level data
#' @rdname ParentOffspring-class
#' @export
ParentOffspring <- function(id=character(),
                            father=character(),
                            mother=character(),
                            offspring=character(),
                            parsedPath=character()){
  new("ParentOffspring", id=id,
      father=father,
      mother=mother,
      offspring=offspring,
      parsedPath=parsedPath)
}

setValidity("ParentOffspring", function(object){
  msg <- TRUE
  if(!all(file.exists(parsedPath(object)))){
    msg <- "Not all source files exist. See parsedPath(object)."
  }
  msg
})

#' @param object a \code{ParentOffspring} object
#' @aliases pedigreeName,ParentOffspring-method
#' @rdname ParentOffspring-class
setMethod("pedigreeName", "ParentOffspring", function(object) object@id)

#' @aliases father,ParentOffspring-method
#' @rdname ParentOffspring-class
setMethod("father", "ParentOffspring", function(object) object@father)

#' @aliases mother,ParentOffspring-method
#' @rdname ParentOffspring-class
setMethod("mother", "ParentOffspring", function(object) object@mother)

#' @aliases offspring,ParentOffspring-method
#' @rdname ParentOffspring-class
setMethod("offspring", "ParentOffspring", function(object) object@offspring)

setMethod("parsedPath", "ParentOffspring", function(object) object@parsedPath)

setMethod("show", "ParentOffspring", function(object){
  cat("Pedigree ID:", pedigreeName(object), "\n")
  cat("father     :", father(object), "\n")
  cat("mother     :", mother(object), "\n")
  cat("offspring  :", paste(offspring(object), collapse=", "), "\n")
})



#' Constructor for ParentOffspringList class
#'
#' @param pedigrees a list of \code{ParentOffspring} objects
#' @rdname ParentOffspringList-class
#' @export
ParentOffspringList <- function(pedigrees=list()){
  ids <- as.character(sapply(pedigrees, pedigreeName))
  new("ParentOffspringList", id=ids, pedigrees=pedigrees)
}

#' @param object a \code{ParentOffspringList} object
#' @aliases pedigreeName,ParentOffspringList-method
#' @rdname ParentOffspringList-class
setMethod("pedigreeName", "ParentOffspringList", function(object) object@id)

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

#' @aliases "[",ParentOffspringList,ANY-method
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

#' @param x a \code{ParentOffspring} object
#' @aliases names,ParentOffspring-method
#' @rdname ParentOffspring-class
setMethod("names", "ParentOffspring", function(x){
  c(father(x), mother(x), offspring(x))
})

setMethod("fileName", "ParentOffspringList", function(object, label){
  datadir <- sapply(object@pedigrees, function(x) dirname(parsedPath(x)[1]), USE.NAMES=FALSE)
  file.path(datadir, paste0(pedigreeName(object), "_", label, ".rds"))
})

setMethod("fileName", "ParentOffspring", function(object, label){
  dirs <- dirname(parsedPath(object))
  ids <- paste0(names(object), "_", label, ".rds")
  file.path(dirs, ids)
})
setMethod("lrrFile", "ParentOffspring", function(object, label="lrr"){
  fileName(object, label)
})

setMethod("bafFile", "ParentOffspring", function(object, label="baf"){
  fileName(object, label)
})

setMethod("gtFile", "ParentOffspring", function(object, label="gt"){
  fileName(object, label)
})



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
