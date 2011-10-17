setMethod("initialize", "LogRratioSet",
	  function(.Object,
		   logRRatio=new("matrix"),
		   BAF=matrix(NA, nrow(logRRatio), ncol(logRRatio)),
		   ...){
		  callNextMethod(.Object,
				 logRRatio=logRRatio,
				 BAF=BAF, ...)
	  })

##setAs("BeadStudioSet", "LogRratioSet",
##      function(from, to){
##	      new("LogRratioSet",
##		  logRRatio=lrr(from),
##		  BAF=baf(from),
##		  phenoData=phenoData(from),
##		  featureData=featureData(from),
##		  protocolData=protocolData(from),
##		  annotation=annotation(from))
##      })



setMethod("logR", "LogRratioSet", function(object) assayDataElement(object, "logRRatio"))
setMethod("lrr", "LogRratioSet", function(object) assayDataElement(object, "logRRatio"))

setReplaceMethod("logR", c("LogRratioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "logRRatio", value)
	 })
setReplaceMethod("lrr", c("LogRratioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "logRRatio", value)
	 })

##setReplaceMethod("logR", c("LogRratioSet", "ffdf"),
##		 function(object, value) {
##			 assayDataElementReplace(object, "logRRatio", value)
##	 })

setMethod("baf", "LogRratioSet",
	  function(object) {
		  assayDataElement(object, "BAF")
	 })

##setReplaceMethod("baf", c("LogRratioSet", "ffdf"),
##		 function(object, value) {
##			 assayDataElementReplace(object, "BAF", value)
##	 })
##setReplaceMethod("baf", c("LogRratioSet", "matrix"),
##		 function(object, value) {
##			 assayDataElementReplace(object, "BAF", value)
##	 })
setReplaceMethod("baf", c("LogRratioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "BAF", value)
	 })
