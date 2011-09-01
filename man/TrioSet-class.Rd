\name{TrioSet-class}
\Rdversion{1.1}
\docType{class}
\alias{TrioSet-class}
\alias{[,TrioSet-method}
\alias{baf<-,TrioSet-method}
\alias{dim,TrioSet-method}
\alias{fatherNames,TrioSet-method}
\alias{logR<-,TrioSet-method}
\alias{mad,TrioSet-method}
\alias{mad<-,TrioSet-method}
\alias{mindist,TrioSet-method}
\alias{mindist<-,TrioSet-method}
\alias{motherNames,TrioSet-method}
\alias{offspringNames,TrioSet-method}
\alias{prune,TrioSet,RangedDataCNV-method}
\alias{todf,TrioSet,RangedData-method}
\alias{updateObject,TrioSet-method}
\alias{xsegment,TrioSet-method}
\alias{xyplot,formula,TrioSet-method}

\title{Class \code{"TrioSet"}}
\description{

A TrioSet is a container for storing high throughput data from
genotyping arrays for case-parent trio designs.  The assay data
elements are 3-d arrays: marker x trio x trio membership.

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TrioSet", assayData, phenoData, featureData, experimentData, annotation, protocolData, ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{phenoData2}:}{Object of class \code{"array"} ~~ }
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
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
No methods defined with class "TrioSet" in the signature.
}
\references{
}
\author{
R. Scharpf
}
\note{

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
