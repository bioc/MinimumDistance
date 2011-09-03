\name{TrioSet-class}
\Rdversion{1.1}
\docType{class}
\alias{TrioSet-class}
\alias{baf<-,TrioSet-method}
\alias{[,TrioSet-method}
\alias{calculateMindist,TrioSet-method}
\alias{dim,TrioSet-method}
\alias{fatherNames,TrioSet-method}
\alias{logR<-,TrioSet-method}
\alias{mad<-,TrioSet-method}
\alias{mad,TrioSet-method}
\alias{mindist<-,TrioSet-method}
\alias{mindist,TrioSet-method}
\alias{motherNames,TrioSet-method}
\alias{offspringNames,TrioSet-method}
\alias{phenoData2<-,TrioSet-method}
\alias{phenoData2,TrioSet-method}
\alias{prune,TrioSet,RangedDataCNV-method}
\alias{todf,TrioSet,RangedData-method}
\alias{updateObject,TrioSet-method}
\alias{xsegment,TrioSet-method}
\alias{xyplot,formula,TrioSet-method}

\title{Class \code{"TrioSet"}}
\description{

  A TrioSet is a container for storing high throughput data assay data
  and metadata from genotyping arrays for case-parent trio designs.  It
  extends eSet.  The primary difference for our application is that the
  assay data elements are 3-d arrays: marker x trio x individual.
  Typically, a TrioSet instance stores the data for a single chromosome
  and an instance of TrioSetList is a list of TrioSets.

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TrioSet", phenoData2, ...)}.

}
\section{Slots}{
  \describe{
    \item{\code{phenoData2}:}{Object of class \code{"arrayOrNULL"} ~~ }
    \item{\code{mindist}:}{Object of class \code{"matrixOrNULL"} ~~ }
    \item{\code{mad}:}{Object of class \code{"matrix"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{LogRatioSet}, directly.
Class \code{"\linkS4class{eSet}"}, by class "LogRatioSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "LogRatioSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "LogRatioSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{baf<-}{\code{signature(object = "TrioSet")}: ... }
    \item{[}{\code{signature(x = "TrioSet")}: ... }
    \item{calculateMindist}{\code{signature(object = "TrioSet")}: ... }
    \item{dim}{\code{signature(x = "TrioSet")}: ... }
    \item{fatherNames}{\code{signature(object = "TrioSet")}: ... }
    \item{logR<-}{\code{signature(object = "TrioSet")}: ... }
    \item{mad<-}{\code{signature(object = "TrioSet")}: ... }
    \item{mad}{\code{signature(x = "TrioSet")}: ... }
    \item{mindist<-}{\code{signature(object = "TrioSet")}: ... }
    \item{mindist}{\code{signature(object = "TrioSet")}: ... }
    \item{motherNames}{\code{signature(object = "TrioSet")}: ... }
    \item{offspringNames}{\code{signature(object = "TrioSet")}: ... }
    \item{phenoData2<-}{\code{signature(object = "TrioSet")}: ... }
    \item{phenoData2}{\code{signature(object = "TrioSet")}: ... }
    \item{prune}{\code{signature(object = "TrioSet", ranges = "RangedDataCNV")}: ... }
    \item{todf}{\code{signature(object = "TrioSet", range = "RangedData")}: ... }
    \item{updateObject}{\code{signature(object = "TrioSet")}: ... }
    \item{xsegment}{\code{signature(object = "TrioSet")}: ... }
    \item{xyplot}{\code{signature(x = "formula", data = "TrioSet")}: ... }
	 }
}
\author{
R. Scharpf
}
\seealso{
  	\code{\linkS4class{TrioSetList}}
}
\examples{
\dontrun{
       R <- array(NA, dim=c(10, 2, 3))
       dimnames(R) <- list(letters[1:10], c("A", "B"), c("F", "M", "O"))
       B <- array(runif(10*2*3, 0, 1), dim=dim(R), dimnames=dimnames(R))
       x <- new("TrioSet", logRRatio=R, BAF=B)
       x[1:2, ]
       x[1:2,1]
       x[1:2,1,1]
       x <- x[, , 1] ## should get all fathers
       x <- x[, , 2] ## should get all mothers
       x <- x[, , 3] ## should get all offspring

       R <- array(NA, dim=c(10, 2, 3))
       dimnames(R) <- list(letters[1:10], c("A", "B"), c("F", "M", "O"))
       B <- array(runif(10*2*3, 0, 1), dim=dim(R), dimnames=dimnames(R))
       fD <- annotatedDataFrameFrom(R, byrow=T)

       ## One could also use 'ff' objects in the assayData insteadof matrices

       ldPath("~/")
       R <- createFF("r", dim=c(10,2,3), vmode="double")
       B <- createFF("r", dim=c(10,2,3), vmode="double")
       x2 <- new("TrioSet", logRRatio=R, BAF=B,featureData=fD)
}
showClass("TrioSet")
}
\keyword{classes}
