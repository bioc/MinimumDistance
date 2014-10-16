#' @include help.R
NULL

setOldClass("ff_array")
setOldClass("ff_matrix")
## setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("arrayORff_array", c("array", "ff_array"))

##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ Pedigree Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Deprecated class for storing pedigree data
#'
#' @slot trios a \code{data.frame} with colnames 'F', 'M', and 'O'
#' containing sample identifiers for the father (F), mother (M), and
#' offspring (O).
#' @slot trioIndex a \code{data.frame}
#' @rdname Pedigree-class
#' @export
setClass("Pedigree", representation(trios="data.frame",
				    trioIndex="data.frame"))


##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ TrioSet Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Deprecated class for storing low-level genomic data for trios
#'
#' This class is deprecated and will be defunct in a future release.
#'
#' @slot fatherPhenoData \code{AnnotatedDataFrame} containing covariates for the father
#' @slot motherPhenoData \code{AnnotatedDataFrame} containing covariates for the mother
#' @slot pedigree an object of class \code{Pedigree}
#' @slot mindist a numeric matrix of the minimum distance for each trio, or NULL
# @slot assayData
# @slot experimentData
# @slot genome
# @slot phenoData
# @slot protocolData
#' @rdname TrioSet-class
#' @docType class
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

#' Deprecated class for storing low-level genomic data for trios
#'
#' This class is deprecated and will be defunct in a future release.
#'
#' @slot fatherPhenoData \code{AnnotatedDataFrame} containing covariates for the father
#' @slot motherPhenoData \code{AnnotatedDataFrame} containing covariates for the mother
#' @slot pedigree an object of class \code{Pedigree}
#' @rdname TrioSetList-class
#' @export
setClass("TrioSetList", contains="gSetList",
	 representation(pedigree="Pedigree",
			fatherPhenoData="AnnotatedDataFrame",
			motherPhenoData="AnnotatedDataFrame"))

setClass("Pedigree2", contains="DataFrame")

## HmmTrioParam
setClass("PennParam", representation(transitionProb="matrix",
                                     transitionNM="numeric",
                                     initialStateProb="numeric",
                                     table1="numeric",
                                     table3="array",
                                     states="matrix",
                                     names="character",
                                     referenceState="character",
                                     prNonMendelian="numeric"))
                                     ##minimum_distance_threshold="numeric",
                                     ##minimum_MAD="numeric",
                                     ##minimum_emission="numeric"))

## MDParam
##setClass("MDParam", representation(


setClass("DNAcopyParam", representation(alpha="numeric",
                                        min.width="integer",
                                        undo.splits="character",
                                        undo.SD="numeric"))

#' Class and methods for parameters of minimum distance algorithm
#'
#' Contains parameters used for circular binary segmentation (package
#' DNAcopy), parameters in the PennCNV hidden Markov model, and
#' parameters used for computing emission probabilities.
#' @slot nMAD a length-one numeric vector
#' @slot dnacopy an object of class \code{DNAcopyParam}
#' @slot penncnv an object of class \code{PennParam}
#' @slot emission an object of class \code{EmissionParam}
#' @slot thin a length-one non-negative integer
#' @export
#' @rdname MinDistParam-class
setClass("MinDistParam", representation(nMAD="numeric",
                                        dnacopy="DNAcopyParam",
                                        penncnv="PennParam",
                                        emission="EmissionParam",
                                        thin="integer"))

#' Object containing the sample identifiers for members in a pedigree
#'
#' Container for registering sample identifiers with membership in a
#' pedigree.  For representing multiple pedigrees, see
#' \code{\linkS4class{ParentOffspringList}}.
#'
#' @slot id length-one character vector providing a family-level id
#' @slot father length-one character vector providing sample ids for father
#' @slot mother length-one character vector providing sample ids for mother
#' @slot offspring character vector providing sample ids for offspring (can have length greater than one if there is more than one offspring)
#' @slot parsedPath character vector providing path to parsed files of the marker-level summaries
#' @rdname ParentOffspring-class
#' @export
#' @examples
#' ParentOffspring()
#' @seealso ParentOffspringList-class
setClass("ParentOffspring", representation(id="character", ## id for pedigree
                                           father="character",
                                           mother="character",
                                           offspring="character",
                                           parsedPath="character"))

#' A list of \code{ParentOffspring} objects
#'
#' Each element of the list is an element of class
#' \code{\linkS4class{ParentOffspring}}.
#'
#' @slot id a character vector of identifiers for the
#' pedigrees. \code{id} must have the same length as \code{pedigrees}
#' @slot pedigrees  A list of ParentOffspring objects.
#' @examples
#' ParentOffspringList()
#' @export
#' @rdname ParentOffspringList-class
setClass("ParentOffspringList", representation(id="character", pedigrees="list"))

#' A container for storing segmentation data for members in a
#' \code{ParentOffspring} family
#'
#' @slot mindist a \code{GRangesList} object
#' @slot offspring a \code{GRangesList} object
#' @slot father a \code{GRanges} object
#' @slot mother a \code{GRanges} object
#' @slot pedigree a \code{ParentOffspring} object
#' @examples
#' data(md_gr)
#' offspring(md_gr)
#' father(md_gr)
#' mother(md_gr)
#' mindist(md_gr)
#' @rdname MinDistGRanges-class
#' @export
setClass("MinDistGRanges", representation(mindist="GRangesList",
                                          offspring="GRangesList",
                                          father="GRanges",
                                          mother="GRanges",
                                          pedigree="ParentOffspring"))


#' Class and methods for MinDistExperiment
#'
#' @slot mindist a matrix
#' @slot pedigree a \code{ParentOffspring} object
#' @export
#' @rdname MinDistExperiment-class
setClass("MinDistExperiment", contains="SnpArrayExperiment",
         representation(mindist="matrix", pedigree="ParentOffspring"))

##setClass("FileViews", representation(path="character",
##                                     pedigree="list",
##                                     cnvar="character",
##                                     bafvar="character",
##                                     fid="character",
##                                     importfun="function",
##                                     annot_pkg="character"))

#' A \code{GRanges}-derived class
#'
#' Contains maximum a posteriori estimates for each genomic interval
#' @examples
#' MDRanges()
#' @rdname MDRanges-class
#' @export
setClass("MDRanges", contains="GRanges")


setClass("ILimit", contains="IRanges")

#' Container for the segmentation results from a MinDistExperiment
#'
#' MinDistPosterior is a \code{GRangesList}-derived container for the
#' segmentation and maximum a posteriori trio copy number states.
#' @export
#' @seealso \code{\link{denovo}}
#' @rdname MinDistPosterior-class
setClass("MinDistPosterior", representation(granges="GRangesList"))

#' A class for filtering genomic intervals called by MinimumDistance
#'
#' Options for filtering include the number of markers spanned by a
#' segment, the posterior probability of the maximum a posteriori
#' estimate of the trio copy number state, and the trio copy number
#' state.
#'
#' @seealso denovo
#' @export
#' @rdname FilterParamMD-class
setClass("FilterParamMD", contains="FilterParam")
