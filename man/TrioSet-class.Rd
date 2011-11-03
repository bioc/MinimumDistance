\name{TrioSet-class}
\Rdversion{1.1}
\docType{class}
\alias{TrioSet-class}
\alias{baf<-,TrioSet-method}
\alias{[,TrioSet-method}
\alias{calculateMindist,TrioSet-method}
\alias{computeBayesFactor,TrioSet-method}
\alias{dim,TrioSet-method}
\alias{fatherNames<-,TrioSet,character-method}
\alias{fatherNames,TrioSet-method}
\alias{initialize,TrioSet-method}
\alias{mad.marker<-,TrioSet,ANY-method}
\alias{mad.mindist<-,TrioSet-method}
\alias{mad.sample<-,TrioSet,ANY-method}
\alias{mindist<-,TrioSet-method}
\alias{mindist,TrioSet-method}
\alias{motherNames<-,TrioSet,character-method}
\alias{motherNames,TrioSet-method}
\alias{offspringNames<-,TrioSet,character-method}
\alias{offspringNames,TrioSet-method}
\alias{phenoData2,TrioSet-method}
\alias{prune,TrioSet,RangedDataCNV-method}
\alias{updateObject,TrioSet-method}
%\alias{xsegment,TrioSet,Pedigree-method}
\alias{xyplot,formula,TrioSet-method}

\title{Class \code{"TrioSet"}}
\description{

    A TrioSet is a container for storing high throughput assay data and
  metadata from genotyping arrays when the study design is case-parent
  trios. In our application, de novo copy number alterations in affected
  offspring were of the primary interest.  Examination of the joint
  distribution of the log R ratios and B allele frequencies across
  members in the trios was a motivating factor in developing the TrioSet
  class. In particular, the assay data elements for a \code{TrioSet}
  object are 3-dimensional arrays: marker x trio x individual.
  Typically, a TrioSet instance stores the data for a single chromosome
  and an instance of TrioSetList is a list of TrioSets.  While having a
  single TrioSet for the entire dataset would simplify the classes of
  objects defined in this package, multiple arrays with thousands of
  trios and roughly a million markers are impractical on most machines.

}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TrioSet", logRRatio, BAF, phenoArray, mindist, mad, ...)}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
    \item{\code{phenoData2}:}{Object of class \code{"array"} ~~ }
    \item{\code{mindist}:}{Object of class \code{"matrix"} ~~ }
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
Class \code{"LogRratioSet"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "LogRratioSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "LogRratioSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "LogRratioSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{baf<-}{\code{signature(object = "TrioSet")}: accessor for B allele frequencies }
    \item{[}{\code{signature(x = "TrioSet")}: ...}
    \item{calculateMindist}{\code{signature(object = "TrioSet")}: Calculate the minimum distance.}
    \item{computeBayesFactor}{\code{signature(object = "TrioSet")}: Compute the maximum a posterior probabilities for the trio copy number states }
    \item{dim}{\code{signature(x = "TrioSet")}: Returns the dimension of the assay data elements. Each assay data element is a three dimensional array with dimensions for features, trios, and family membership, respectively.}
    \item{fatherNames<-}{\code{signature(object = "TrioSet", value = "character")}: Assign father names}
    \item{fatherNames}{\code{signature(object = "TrioSet")}: }
    \item{initialize}{\code{signature(.Object = "TrioSet")}: ... }
    \item{mad.marker<-}{\code{signature(x = "TrioSet", value = "ANY")}: ...}
    \item{mad.mindist<-}{\code{signature(x = "TrioSet")}: ... }
    \item{mad.sample<-}{\code{signature(x = "TrioSet", value = "ANY")}: ... }
    \item{mindist<-}{\code{signature(object = "TrioSet")}: Assign the minimum distance to slot \code{mindist} }
    \item{mindist}{\code{signature(object = "TrioSet")}: Accessor for the minimum distance}
    \item{motherNames<-}{\code{signature(object = "TrioSet", value = "character")}: ... }
    \item{motherNames}{\code{signature(object = "TrioSet")}: ... }
    \item{offspringNames<-}{\code{signature(object = "TrioSet", value = "character")}: ... }
    \item{offspringNames}{\code{signature(object = "TrioSet")}: ... }
    \item{phenoData2}{\code{signature(object = "TrioSet")}:  accessor
      for phenoData2 slot.  }
    \item{prune}{\code{signature(object = "TrioSet", ranges = "RangedDataCNV")}: ... }
    \item{xsegment}{\code{signature(object = "TrioSet", pedigreeData = "Pedigree")}: Run circular binary segmentation on the minimum distance}
    \item{xyplot}{\code{signature(x = "formula", data = "TrioSet")}: plot the log R ratios and B allele frequencies for a trio}
  }
}
\author{R. Scharpf}


\seealso{
	\code{\linkS4class{TrioSetList}}
}
\examples{
	showClass("TrioSet")
	## instantiate a TrioSet with an array of log R ratios (logRR) and B
	## allele frequencies
	logRR <- array(rnorm(10*2*3), dim=c(10, 2, 3))
	dimnames(logRR) <- list(paste("marker", 1:10, sep=""), paste("trio", 1:2, sep=""), c("F", "M", "O"))
	BAFs <- array(runif(10*2*3, 0, 1), dim=dim(logRR), dimnames=dimnames(logRR))
	trioSet <- new("TrioSet", logRRatio=logRR, BAF=BAFs)
	## log R ratios for first five markers for all trios
	lrr(trioSet)[1:5, , ]
	## log R ratios for first five markers and the first trio
	lrr(trioSet)[1:5, 1, ]
	## log R ratios for only the offspring in the first trio
	lrr(trioSet)[1:5, 1, "O"]
	## similarly for B allele frequencies
	baf(trioSet)[1:5, 1, "O"]
	if(require("ff")){
		## for large datasets, one can generate ff arrays for the assay data
		## rather than store these in memory
		##
		## First, indicate the directory to store the ff objects
		\dontrun{
			ldPath(".")
		}
		logRR.ff <- oligoClasses:::initializeBigArray(name="lrr", dim=c(10,2,3), vmode="double")
		## after initializion, one can read subsets of the data and write
		## to disk by the "[<-" assignment method
		## For example, here we write the data for 10 markers to disk
		logRR.ff[1:10, , ] <- logRR[1:10, , ]
		## Similarly for BAFs
		BAF.ff <- oligoClasses:::initializeBigArray(name="baf", dim=c(10,2,3), vmode="double")
		BAF.ff[1:10, , ] <- BAFs[1:10, , ]
		trioSet.ff <- new("TrioSet", logRRatio=logRR.ff, BAF=BAF.ff)
	}
}
\keyword{classes}
