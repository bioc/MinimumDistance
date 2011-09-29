setMethod("offspringNames", signature(object="Pedigree"), function(object) object$O)
setMethod("fatherNames", signature(object="Pedigree"), function(object) object$F)
setMethod("motherNames", signature(object="Pedigree"), function(object) object$M)
