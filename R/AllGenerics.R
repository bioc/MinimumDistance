#' @include help.R
NULL

#' Getter and setter for the minimum distance statistic
#'
#' @param object see \code{showMethods("mindist")}
#' @rdname mindist
#' @export
setGeneric("mindist", function(object) standardGeneric("mindist"))

#' @param value a matrix of the minimum distance
#' @rdname mindist
#' @export
setGeneric("mindist<-", function(object,value) standardGeneric("mindist<-"))

setGeneric("range.index", function(object) standardGeneric("range.index"))

#' Deprecated function to calculate the minimum distance
#'
#'   The 'minimum distance' is the minimum signed absolute difference of
#'    the parental log R ratios and the offspring log R
#'    ratios. Specifically, let |O-F| denote the absolute difference in
#'    the log R ratios comparing offspring to father and |O-M| the
#'    absolute difference in the log R ratios comparing offspring to
#'    mother.  The minimum distance at a marker is the signed minimum of
#'    |O-M| and |O-F|.  After segmentation of the minimum distance,
#'    non-zero segments can indicate a de novo difference in the log R
#'    ratio of the offspring and either parent.  For example, a positive
#'    minimum distance suggests that the log R ratio from the offspring is
#'    greater than the log R ratio of either parent.
#' @param object see \code{showMethods("calculateMindist")}
#' @param ... Ignored
#' @export
setGeneric("calculateMindist", function(object, ...) standardGeneric("calculateMindist"))

setGeneric("mad")
setGeneric("stack")

setGeneric("prune", function(object,  ranges, ...) standardGeneric("prune"))
setGeneric("computeBayesFactor", function(object, ranges, ...) standardGeneric("computeBayesFactor"))
setGeneric("todf", function(object, rangeData, frame, ...) standardGeneric("todf"))

#' Deprecated functions and methods
#'
#' These functions will be defunct in a future release.
#'
#' @param object see \code{showMethods("offspringNames")}
#' @rdname Deprecated
#' @export
setGeneric("offspringNames", function(object) standardGeneric("offspringNames"))

#' @param value a character vector of offspring identifiers
#' @rdname Deprecated
#' @export
setGeneric("offspringNames<-", function(object,value) standardGeneric("offspringNames<-"))
setGeneric("fatherNames", function(object) standardGeneric("fatherNames"))
setGeneric("fatherNames<-", function(object,value) standardGeneric("fatherNames<-"))
setGeneric("motherNames", function(object) standardGeneric("motherNames"))
setGeneric("motherNames<-", function(object,value) standardGeneric("motherNames<-"))

setGeneric("allNames", function(object) standardGeneric("allNames"))

#' @export
#' @rdname Deprecated
setGeneric("trios", function(object) standardGeneric("trios"))

setGeneric("minimumDistance", function(object, ...) standardGeneric("minimumDistance"))
setGeneric("trioplot", function(formula, object, range, ...) standardGeneric("trioplot"))

#' A wrapper for DNAcopy's segment function
#'
#' Methods for circular binary segmentation.
#'
#' @param object see \code{showMethods{segment2}}
#' @param ... Additional arguments passed to DNAcopy's \code{segment}.
#' @export
#' @seealso \code{\link{segment}}
#' @rdname segment2
setGeneric("segment2", function(object, ...) standardGeneric("segment2"))

#' Deprecated wrapper for computing the median absolute deviation of
#' low-level summaries
#'
#' @param object see \code{showMethods("mad2")}
#' @param byrow logical if TRUE, compute the median absolute deviation of the rows of a matrix
#' @param ... additional arguments to \code{\link{mad}}
#' @export
setGeneric("mad2", function(object, byrow=FALSE, ...) standardGeneric("mad2"))

setGeneric("assayDataList", function(object) standardGeneric("assayDataList"))
setGeneric("fatherPhenoData", function(object) standardGeneric("fatherPhenoData"))
setGeneric("motherPhenoData", function(object) standardGeneric("motherPhenoData"))
setGeneric("offspringPhenoData", function(object) standardGeneric("offspringPhenoData"))
setGeneric("featureDataList", function(object) standardGeneric("featureDataList"))
setGeneric("trioIndex", function(object) standardGeneric("trioIndex"))
setGeneric("chromosomeList", function(object) standardGeneric("chromosomeList"))
setGeneric("sampleNames2", function(object) standardGeneric("sampleNames2"))

## defunct
setGeneric("gcSubtract", function(object, ...) setGeneric("gcSubtract"))

#' Computes the maximum a posteriori trio copy number state for the
#' segmented minimum distance
#'
#' This functions is deprecated and will be defunct in a future
#' release.  The replacement function is MAP2.
#'
#' @param object see \code{showMethods(MAP)}
#' @param ranges A \code{GRanges} object
#' @param id character string for sample identifier
#' @param TAUP scalar for transition probabilities
#' @param tauMAX the maximum probability that the current state is the same as the previous state
#' @param cnStates character vector for hidden Markov model state labels
#' @param pr.nonmendelian numeric:  the a priori probability of a non-Mendelian copy number alteration
#' @param mdThr a length-one numeric vector.  A minimum distance below this threshold in absolute value will not be evaluated for copy number alterations.
#' @param ...  Ignored.
#' @rdname MAP
#' @export
setGeneric("MAP", function(object, ranges, id,
			   TAUP=1e10, tauMAX=1-5e-8,
			   cnStates=c(-2, -0.4, 0, 0, 0.4, 1),
			   pr.nonmendelian=1.5e-6, mdThr=0.9, ...) standardGeneric("MAP"))

setGeneric("table1", function(object) standardGeneric("table1"))
setGeneric("table3", function(object) standardGeneric("table3"))
setGeneric("stateNames", function(object) standardGeneric("stateNames"))
setGeneric("referenceState", function(object) standardGeneric("referenceState"))
setGeneric("prNonMendelian", function(object) standardGeneric("prNonMendelian"))
setGeneric("minimum_distance_threshold", function(object) standardGeneric("minimum_distance_threshold"))

setGeneric("initialStateProb", function(object) standardGeneric("initialStateProb"))
setGeneric("initialStateProb<-", function(object,value) standardGeneric("initialStateProb<-"))
setGeneric("transitionProb", function(object) standardGeneric("transitionProb"))

setGeneric("minimum_MAD", function(object) standardGeneric("minimum_MAD"))

#' Setter and getter for number of median absolute deviations the mean
#' minimum distance of a genomic interval is from zero
#'
#' @param object see \code{showMethods("nMAD")}
#' @rdname nMAD
#' @export
setGeneric("nMAD", function(object) standardGeneric("nMAD"))

#' @param value a length-one numeric vector
#' @rdname nMAD
#' @export
setGeneric("nMAD<-", function(object,value) standardGeneric("nMAD<-"))

## export
## setGeneric("xyplot") ##, signature=c("x", "data"))

setGeneric("alpha", function(object) standardGeneric("alpha"))
setGeneric("min.width", function(object) standardGeneric("min.width"))
setGeneric("undo.splits", function(object) standardGeneric("undo.splits"))
setGeneric("undo.SD", function(object) standardGeneric("undo.SD"))

setGeneric("subsetAndSort", function(object, autosomes) standardGeneric("subsetAndSort"))

#' @rdname ParentOffspring-class
#' @export
setGeneric("offspring", function(object) standardGeneric("offspring"))

#' @rdname ParentOffspring-class
#' @export
setGeneric("mother", function(object) standardGeneric("mother"))

#' Sample identifier for the father in a pedigree
#'
#' Accessor for the sample identifiers for the members in a
#' pedigree
#' @rdname ParentOffspring-class
#' @export
setGeneric("father", function(object) standardGeneric("father"))

setGeneric("acfs", function(x) standardGeneric("acfs"))
setGeneric("penncnv", function(object) standardGeneric("penncnv"))

#' @rdname Deprecated
#' @export
setGeneric("pedigree", function(object) standardGeneric("pedigree"))

#' @rdname Deprecated
#' @export
setGeneric("pedigree<-", function(object,value) standardGeneric("pedigree<-"))

#' Accessor for pedigree name
#'
#' @param object a \code{ParentOffspring} or \code{ParentOffspringList} object
#' @export
#' @seealso \code{\linkS4class{ParentOffspring}} \code{\linkS4class{ParentOffspringList}}
setGeneric("pedigreeName", function(object) standardGeneric("pedigreeName"))

#' Computes maximum a posteriori estimate for the trio copy number state
#'
#' @param object An object of class \code{MinDistExperiment}
#' @param mdgr An object of class \code{MinDistGRanges}, \code{GRangesList}, or \code{GRanges}.
#' @param param An object of class \code{MinDistParam}.
#' @param ... ignored
#' @return An object of class \code{MinDistPosterior}
#' @examples
#'   library(oligoClasses)
#'   library(VanillaICE)
#'   ## A MinDistExperiment object:
#'   data(md_exp)
#'   ## Segmented data
#'   data(md_gr)
#'   e_param <- EmissionParam(temper=1, p_outlier=1/100)
#'   param <- MinDistParam(thin=1L, emission=e_param)
#' \dontrun{
#'   md_g <- MAP2(md_exp, md_gr, param)
#' }
#' @export
setGeneric("MAP2", function(object, mdgr, param=MinDistParam(), ...) standardGeneric("MAP2"))

setGeneric("posteriorLogOdds", function(object) standardGeneric("posteriorLogOdds"))
setGeneric("posteriorLogRR", function(object) standardGeneric("posteriorLogRR"))
setGeneric("thin", function(object) standardGeneric("thin"))

#' Constructor for \code{MinDistExperiment} class
#'
#' @param object see \code{showMethods(MinDistExperiment)}
#' @param pedigree a \code{ParentOffspring} object
#' @param ... ignored
#' @return an object of class \code{MinDistExperiment}
#' @export
setGeneric("MinDistExperiment", function(object=ArrayViews(), pedigree=ParentOffspring(), ...) standardGeneric("MinDistExperiment"))

setGeneric("is221", function(object) standardGeneric("is221"))


setGeneric("mdgrFile", function(object, label="mdgr") standardGeneric("mdgrFile"))
setGeneric("mapFile", function(object, label="mdgr") standardGeneric("mapFile"))
setGeneric("seq_along2", function(along.with) standardGeneric("seq_along2"))

setGeneric("computePosterior", function(object, granges, param) standardGeneric("computePosterior"))


#' Methods for filtering MinDistExperiment objects
#'
#' Filter a MinDistExperiment object to exclude markers with missing
#' values in the low-level summaries, exclude markers that lie in
#' segments (\code{granges} argument) with small minimum distance
#' values (unlikely to be de novo)
#'
#' @return a \code{MinDistExperiment}
#' @param object A \code{MinDistExperiment}
#' @param granges A \code{GRanges}, \code{GRangesList}, or \code{MinDistGRanges} object
#' @param param a \code{MinDistParam} object
#' @export
setGeneric("filterExperiment", function(object, granges, param) standardGeneric("filterExperiment"))

#' @rdname denovo
#' @export
setGeneric("denovoHemizygous", function(object, filters=FilterParamMD(state="221")) standardGeneric("denovoHemizygous"))


#' @rdname denovo
#' @export
setGeneric("denovoHomozygous", function(object, filters=FilterParamMD(state="220")) standardGeneric("denovoHomozygous"))

#' Filter the genomic intervals for denovo copy number states
#'
#' This function filters the genomic intervals for denovo events.
#'
#'  The function \code{denovo} filters genomic intervals for states
#' '220', '221', and '224', corresponding to denovo homozygous
#' deletion, denovo hemizygous deletion, and denovo duplication,
#' respectively.
#'
#' \code{denovoHemizygous} filters genomic intervals for state '221'.
#'
#' \code{denovoHomozygous} filters genomic intervals for state '220'
#' @param object see \code{showMethods(denovo)} for a list of defined methods
#' @param filters an object of class \code{FilterParamMD}
#' @rdname denovo
#' @export
#' @seealso FilterParamMD-class
setGeneric("denovo", function(object, filters=FilterParamMD(state=c("220", "221", "224"))) standardGeneric("denovo"))

#' @rdname denovo
#' @export
setGeneric("denovoDuplication", function(object, filters=FilterParamMD(state="224")) standardGeneric("denovoDuplication"))

#' Plot marker-level summaries for a genomic interval of interest
#'
#' @param object see \code{showMethods("plotDenovo")}
#' @param g a \code{MDRanges} object
#' @param param a \code{HmmTrellisParam} object
#' @export
#' @rdname plotDenovo
setGeneric("plotDenovo",
           function(object, g, param) standardGeneric("plotDenovo"))
