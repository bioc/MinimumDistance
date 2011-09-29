setMethod("sampleSheet", signature(object="TrioAnnotation"),
	  function(object) object@sampleSheet)
setMethod("pedigree", signature(object="TrioAnnotation"),
	  function(object) object@pedigree)
setMethod("sampleNames", signature(object="TrioAnnotation"),
	  function(object) sampleNames(sampleSheet(object)))
setMethod("nrow", signature(x="TrioAnnotation"),
	  function(x) nrow(pedigree(x)))
setMethod("offspringNames", signature(object="TrioAnnotation"), function(object){
	offspringNames(pedigree(object))
})
setMethod("fatherNames", signature(object="TrioAnnotation"), function(object){
	fatherNames(pedigree(object))
})
setMethod("motherNames", signature(object="TrioAnnotation"), function(object){
	motherNames(pedigree(object))
})

setAs("TrioAnnotation", "array",
      function(from){
	      object <- from
	      sample.sheet <- sampleSheet(object)
	      ss <- array(NA, dim=c(nrow(object), ncol(sample.sheet), 3),
			  dimnames=list(offspringNames(object),
			  colnames(sample.sheet),
			  colnames(pedigree(object))))
	      father.index <- match(fatherNames(object),
				    sampleNames(object))
	      mother.index <- match(motherNames(object),
				    sampleNames(object))
	      offspring.index <- match(offspringNames(object),
				       sampleNames(object))
	      ss[, , "F"] <- as.matrix(sample.sheet[father.index, ])
	      ss[, , "M"] <- as.matrix(sample.sheet[mother.index, ])
	      ss[, , "O"] <- as.matrix(sample.sheet[offspring.index, ])
	      return(ss)
	  })

setMethod("initialize", signature(.Object="TrioAnnotation"),
	  function(.Object,
		   pedigree=new("Pedigree"),
		   sampleSheet=new("SampleSheet")){
		  .Object@sampleSheet <- sampleSheet
		  .Object@pedigree <- pedigree
		  return(.Object)
	  })
