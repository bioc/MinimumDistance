setOldClass("ff_array")
setOldClass("ff_matrix")
##setClass("LogRratioSet", contains="eSet")
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("arrayORff_array", c("array", "ff_array"))
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ Pedigree Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setClass("Pedigree", representation(trios="data.frame",
				    trioIndex="data.frame"))


##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ TrioSet Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setClass("TrioSet", contains="gSet",
	 representation(fatherPhenoData="AnnotatedDataFrame",
			motherPhenoData="AnnotatedDataFrame",
			pedigree="Pedigree",
			mindist="matrixOrNULL"))
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ TrioSetList Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setClass("TrioSetList", contains="gSetList",
	 representation(pedigree="Pedigree",
			fatherPhenoData="AnnotatedDataFrame",
			motherPhenoData="AnnotatedDataFrame"))

setClass("Pedigree2", contains="DataFrame")


#' @importClassesFrom VanillaICE SnpArrayExperiment SnpGRanges
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
