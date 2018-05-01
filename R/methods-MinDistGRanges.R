#' Constructor for \code{MinDistGRanges} class
#'
#' The \code{MinDistGRanges} class contains the segmentation of the
#' father, mother, offspring, and the minimum distance for each
#' possible parent-offspring trio.  For the parents, the segmentation
#' results are expected to be in \code{GRanges} format.  To accomodate
#' multiple-offspring families, both the offspring segments and
#' minimum distance segments should be of class \code{GRangesList}
#' where the length of the list corresponds to the number of
#' offspring.
#'
#' @param mindist a \code{GRangesList} object
#' @param offspring  a \code{GRangesList} object
#' @param father a \code{GRanges} object
#' @param mother a \code{GRanges} object
#' @param pedigree a \code{ParentOffspring} object
#' @examples
#' MinDistGRanges()
#' @export
MinDistGRanges <- function(mindist=GRangesList(),
                           offspring=GRangesList(),
                           father=GRanges(),
                           mother=GRanges(),
                           pedigree=ParentOffspring()){
                           ##mad=numeric(), acf=numeric()){
  new("MinDistGRanges", mindist=mindist, offspring=offspring,
      father=father, mother=mother, pedigree=pedigree)
}

#' @param x a \code{MinDistGRanges} object
#' @aliases names,MinDistGRanges-method
#' @rdname MinDistGRanges-class
setMethod("names", "MinDistGRanges", function(x) x@id)

#' @param object a \code{MinDistGRanges} object
#' @aliases mindist,MinDistGRanges-method
#' @rdname MinDistGRanges-class
setMethod("mindist", "MinDistGRanges", function(object) object@mindist)

#' @param value a \code{GRangesList} object
#' @aliases mindist<-,MinDistGRanges,GRangesList-method
#' @rdname MinDistGRanges-class
setReplaceMethod("mindist", c("MinDistGRanges", "GRangesList"), function(object, value) {
  object@mindist <- value
  object
})

#' @aliases offspring,MinDistGRanges-method
#' @rdname MinDistGRanges-class
setMethod("offspring", "MinDistGRanges", function(object) object@offspring)

#' @aliases mother,MinDistGRanges-method
#' @rdname MinDistGRanges-class
setMethod("mother", "MinDistGRanges", function(object) object@mother)

#' @aliases father,MinDistGRanges-method
#' @rdname MinDistGRanges-class
setMethod("father", "MinDistGRanges", function(object) object@father)


##setMethod("mad", "MinDistGRanges", function(x, center = median(x), constant = 1.4826, na.rm = FALSE, low = FALSE, high = FALSE) x@mad)
##setMethod("acfs", "MinDistGRanges", function(x) x@acf)
##pedigree_id <- function(object) object@pedigree_id

#' @aliases pedigree,MinDistGRanges-method
#' @rdname MinDistGRanges-class
setMethod("pedigree", "MinDistGRanges", function(object) object@pedigree)

#' @aliases show,MinDistGRanges-method
#' @rdname MinDistGRanges-class
setMethod("show", "MinDistGRanges", function(object){
  cat("An object of class 'MinDistGRanges' \n")
  ##cat("  pedigree id:  ", pedigree_id(object), "\n")
  offspr <- offspring(pedigree(object))
  mdr <- mindist(object)
  L <- elementNROWS(offspring(object))
  L2 <- elementNROWS(mdr)
  mdrnames <- names(mdr)
  cat("  mindist:  GRangesList of length", length(offspring(object)), "\n")
  for(i in seq_along(offspr)){
    cat("     o  ", paste0(mdrnames[i], ": ", L2[i], " ranges\n"))
  }
  cat("  offspring:  GRangesList of length", length(offspring(object)), "\n")
  for(i in seq_along(L2)){
    cat("     o  ", paste0(offspr[i], ": ", L[i], " ranges\n"))
  }
  cat("  mother:  GRanges of length", length(mother(object)), "\n")
  cat("  father:  GRanges of length", length(father(object)), "\n")
  show(pedigree(object))
##  mads <- round(mad(object), 2)
##  autocorrelations <- paste(round(acfs(object), 2), collapse=", ")
##  cat("  MADs: ", mads, "\n")
##  cat("  lag10 ACFs: ", autocorrelations, "\n")
})


#' @aliases offspring,GRangesList-method
#' @rdname MinDistGRanges-class
setMethod("offspring", "GRangesList", function(object) {
  object <- object[grep("offspring", names(object))]
})
