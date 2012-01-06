setGeneric("mindist", function(object) standardGeneric("mindist"))
setGeneric("mindist<-", function(object,value) standardGeneric("mindist<-"))
setGeneric("loglik", function(object) standardGeneric("loglik"))
setGeneric("loglik<-", function(object, value) standardGeneric("loglik<-"))
setGeneric("range.index", function(object) standardGeneric("range.index"))
setGeneric("calculateMindist", function(object, ...) standardGeneric("calculateMindist"))
setGeneric("lrr<-", function(object, value) standardGeneric("lrr<-"))
setGeneric("baf<-", function(object, value) standardGeneric("baf<-"))
setGeneric("mad")
setGeneric("ncol")
##setGeneric("mad<-", function(x, value) standardGeneric("mad<-"))
##setGeneric("mad.sample", function(x) standardGeneric("mad.sample"))
##setGeneric("mad.marker", function(x) standardGeneric("mad.marker"))
##setGeneric("mad.sample<-", function(x, value) standardGeneric("mad.sample<-"))
##setGeneric("mad.marker<-", function(x, value) standardGeneric("mad.marker<-"))
##setGeneric("mad.mindist<-", function(x, value) standardGeneric("mad.mindist<-"))
##setGeneric("mad.mindist", function(x) standardGeneric("mad.mindist"))
##setGeneric("phenoData2", function(object) standardGeneric("phenoData2"))
##setGeneric("phenoData2<-", function(object,value) standardGeneric("phenoData2<-"))
##setGeneric("varLabels2", function(object) standardGeneric("varLabels2"))
##setGeneric("prune", function(object, ranges, id,
##			     lambda=0.05,
##			     min.change=0.1,
##			     min.coverage=3,
##			     scale.exp=0.02,
##			     verbose, ...)
##	   standardGeneric("prune"))

setGeneric("prune", function(object,  ranges, ...) standardGeneric("prune"))

setGeneric("computeBayesFactor", function(object,
					  ranges,
##					  mad.marker,
##					  mad.sample,
##					  returnEmission=FALSE,
##					  verbose=TRUE,
					  ...)
	   standardGeneric("computeBayesFactor"))
setGeneric("todf", function(object, rangeData, frame, ...) standardGeneric("todf"))
setGeneric("offspringNames", function(object) standardGeneric("offspringNames"))
setGeneric("offspringNames<-", function(object,value) standardGeneric("offspringNames<-"))
setGeneric("fatherNames", function(object) standardGeneric("fatherNames"))
setGeneric("fatherNames<-", function(object,value) standardGeneric("fatherNames<-"))
setGeneric("motherNames", function(object) standardGeneric("motherNames"))
setGeneric("motherNames<-", function(object,value) standardGeneric("motherNames<-"))
##setGeneric("trioNames", function(object) standardGeneric("trioNames"))
##setGeneric("rbind", function(..., deparse.level=1) standardGeneric("rbind"),
##           signature = "...")
setGeneric("width", function(x) standardGeneric("width"))
setGeneric("pedigree", function(object) standardGeneric("pedigree"))
setGeneric("allNames", function(object) standardGeneric("allNames"))
setGeneric("trios", function(object) standardGeneric("trios"))
setGeneric("minimumDistance", function(object, ...) standardGeneric("minimumDistance"))
setGeneric("trioplot", function(formula, object, range, ...) standardGeneric("trioplot"))
setGeneric("segment2", function(object, pos, chrom, id=NULL, featureNames, ...) standardGeneric("segment2"))
setGeneric("mad2", function(object, byrow=FALSE, ...) standardGeneric("mad2"))
##setGeneric("order2", function(object, ...) standardGeneric("order2"))

setGeneric("assayDataList", function(object) standardGeneric("assayDataList"))
setGeneric("fatherPhenoData", function(object) standardGeneric("fatherPhenoData"))
setGeneric("motherPhenoData", function(object) standardGeneric("motherPhenoData"))
setGeneric("offspringPhenoData", function(object) standardGeneric("offspringPhenoData"))
setGeneric("featureDataList", function(object) standardGeneric("featureDataList"))
setGeneric("trios", function(object) standardGeneric("trios"))
setGeneric("trioIndex", function(object) standardGeneric("trioIndex"))
setGeneric("chromosomeList", function(object) standardGeneric("chromosomeList"))
