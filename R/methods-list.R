setMethod("calculateMindist", signature(object="list"),
	  function(object, outdir, ...){
		  foreach(elt=object, .packages=c("ff", "MinimumDistance")) %dopar% calculateMindist(elt, outdir=outdir)
	  })


unstack <- function(object){
	lrrs <- lapply(object, lrr)
	bafs <- lapply(object, baf)
	new("TrioSetList",
	    assayDataList=AssayDataList(BAF=bafs, logRRatio=lrrs),
	    featureDataList=lapply(object, featureData),
	    phenoData=phenoData(object[[1]]),
	    motherPhenoData=motherPhenoData(object[[1]]),
	    fatherPhenoData=fatherPhenoData(object[[1]]),
	    chromosome=sapply(object, function(x) unique(chromosome(x))))
}

