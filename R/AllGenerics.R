#' @include help.R

setGeneric("mindist", function(object) standardGeneric("mindist"))
setGeneric("mindist<-", function(object,value) standardGeneric("mindist<-"))
setGeneric("range.index", function(object) standardGeneric("range.index"))
setGeneric("calculateMindist", function(object, ...) standardGeneric("calculateMindist"))
setGeneric("mad")
setGeneric("ncol")
setGeneric("prune", function(object,  ranges, ...) standardGeneric("prune"))
setGeneric("computeBayesFactor", function(object, ranges, ...) standardGeneric("computeBayesFactor"))
setGeneric("todf", function(object, rangeData, frame, ...) standardGeneric("todf"))
setGeneric("offspringNames", function(object) standardGeneric("offspringNames"))
setGeneric("offspringNames<-", function(object,value) standardGeneric("offspringNames<-"))
setGeneric("fatherNames", function(object) standardGeneric("fatherNames"))
setGeneric("fatherNames<-", function(object,value) standardGeneric("fatherNames<-"))
setGeneric("motherNames", function(object) standardGeneric("motherNames"))
setGeneric("motherNames<-", function(object,value) standardGeneric("motherNames<-"))
setGeneric("pedigree", function(object) standardGeneric("pedigree"))
setGeneric("allNames", function(object) standardGeneric("allNames"))
setGeneric("trios", function(object) standardGeneric("trios"))
setGeneric("minimumDistance", function(object, ...) standardGeneric("minimumDistance"))
setGeneric("trioplot", function(formula, object, range, ...) standardGeneric("trioplot"))
setGeneric("segment2", function(object, ...) standardGeneric("segment2"))
setGeneric("mad2", function(object, byrow=FALSE, ...) standardGeneric("mad2"))
setGeneric("assayDataList", function(object) standardGeneric("assayDataList"))
setGeneric("fatherPhenoData", function(object) standardGeneric("fatherPhenoData"))
setGeneric("motherPhenoData", function(object) standardGeneric("motherPhenoData"))
setGeneric("offspringPhenoData", function(object) standardGeneric("offspringPhenoData"))
setGeneric("featureDataList", function(object) standardGeneric("featureDataList"))
setGeneric("trios", function(object) standardGeneric("trios"))
setGeneric("trioIndex", function(object) standardGeneric("trioIndex"))
setGeneric("chromosomeList", function(object) standardGeneric("chromosomeList"))
setGeneric("sampleNames2", function(object) standardGeneric("sampleNames2"))
setGeneric("gcSubtract", function(object, ...) setGeneric("gcSubtract"))
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
setGeneric("transitionProb", function(object) standardGeneric("transitionProb"))

setGeneric("minimum_MAD", function(object) standardGeneric("minimum_MAD"))
setGeneric("minimum_emission", function(object) standardGeneric("minimum_emission"))

## #' @export
## #' @importFrom lattice xyplot
## setGeneric("xyplot") ##, signature=c("x", "data"))
