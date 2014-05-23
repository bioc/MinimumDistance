#' @include help.R
NULL

setOldClass("ff_array")
setOldClass("ff_matrix")
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("arrayORff_array", c("array", "ff_array"))

##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ Pedigree Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @export
setClass("Pedigree", representation(trios="data.frame",
				    trioIndex="data.frame"))


##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ TrioSet Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @rdname TrioSet-class
#' @export
setClass("TrioSet", contains="gSet",
	 representation(fatherPhenoData="AnnotatedDataFrame",
			motherPhenoData="AnnotatedDataFrame",
			pedigree="Pedigree",
			mindist="matrixOrNULL"))
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ TrioSetList Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @export
setClass("TrioSetList", contains="gSetList",
	 representation(pedigree="Pedigree",
			fatherPhenoData="AnnotatedDataFrame",
			motherPhenoData="AnnotatedDataFrame"))

setClass("Pedigree2", contains="DataFrame")


#' @importClassesFrom VanillaICE SnpArrayExperiment SnpGRanges
#' @export
setClass("TrioExperiment",
         representation(pedigree="Pedigree"),
         contains="SnpArrayExperiment")


##setGeneric("TrioExperiment",
##           function(..., pedigree=Pedigree()))
setMethod("initialize", "TrioExperiment", function(.Object, ..., pedigree=Pedigree()){
  .Object <-  callNextMethod(.Object, ...)
  .Object@pedigree <- pedigree
  .Object
})

TrioExperiment <- function(..., pedigree=Pedigree()){
  new("TrioExperiment", ..., pedigree=pedigree)
}


## HmmTrioParam
setClass("PennParam", representation(transitionProb="matrix",
                                     initialStateProb="numeric",
                                     table1="numeric",
                                     table3="array",
                                     states="matrix",
                                     names="character",
                                     referenceState="character",
                                     prNonMendelian="numeric",
                                     minimum_distance_threshold="numeric",
                                     minimum_MAD="numeric",
                                     minimum_emission="numeric"))

## MDParam
##setClass("MDParam", representation(
