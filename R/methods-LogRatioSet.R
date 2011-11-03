setMethod("initialize", "LogRratioSet",
	  function(.Object,
		   logRRatio=new("matrix"),
		   BAF=matrix(NA, nrow(logRRatio), ncol(logRRatio)),
		   ...){
		  callNextMethod(.Object,
				 logRRatio=logRRatio,
				 BAF=BAF, ...)
	  })
setMethod("lrr", "LogRratioSet", function(object) assayDataElement(object, "logRRatio"))
setReplaceMethod("lrr", c("LogRratioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "logRRatio", value)
	 })
setMethod("baf", "LogRratioSet",
	  function(object) {
		  assayDataElement(object, "BAF")
	 })
setReplaceMethod("baf", c("LogRratioSet", "ANY"),
		 function(object, value) {
			 assayDataElementReplace(object, "BAF", value)
	 })
