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

#' @aliases show,ParentOffspring-method
#' @rdname ParentOffspring-class
setMethod("show", "ParentOffspring", function(object){
  cat("Pedigree ID:", pedigreeName(object), "\n")
  cat("father     :", father(object), "\n")
  cat("mother     :", mother(object), "\n")
  cat("offspring  :", paste(offspring(object), collapse=", "), "\n")
})



#' @param x a \code{ParentOffspring} object
#' @aliases names,ParentOffspring-method
#' @rdname ParentOffspring-class
setMethod("names", "ParentOffspring", function(x){
  c(father(x), mother(x), offspring(x))
})


setMethod("fileName", "ParentOffspring", function(object, label){
  dirs <- dirname(parsedPath(object))
  ids <- paste0(names(object), "_", label, ".rds")
  file.path(dirs, ids)
})

setMethod("lrrFile", "ParentOffspring", function(object){
  fileName(object, "lrr")
})

setMethod("bafFile", "ParentOffspring", function(object){
  fileName(object, "baf")
})

setMethod("gtFile", "ParentOffspring", function(object){
  fileName(object, "gt")
})
