\name{TrioSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{TrioSetList-class}
\alias{allNames,TrioSetList-method}
\alias{[,TrioSetList-method}
\alias{calculateMindist,TrioSetList-method}
\alias{computeBayesFactor,TrioSetList-method}
\alias{fatherNames,TrioSetList-method}
\alias{initialize,TrioSetList-method}
\alias{mad.marker<-,TrioSetList,list-method}
\alias{mad.mindist<-,TrioSetList-method}
\alias{mad.sample<-,TrioSetList,matrix-method}
\alias{mindist<-,TrioSetList-method}
\alias{mindist,TrioSetList-method}
\alias{minimumDistance,TrioSetList-method}
\alias{motherNames,TrioSetList-method}
\alias{ncol,TrioSetList-method}
\alias{nrow,TrioSetList-method}
\alias{offspringNames,TrioSetList-method}
\alias{pedigree,TrioSetList-method}
\alias{prune,TrioSetList,RangedDataCNV-method}
\alias{sampleSheet,TrioSetList-method}
\alias{stack,TrioSetList-method}
\alias{trios,TrioSetList-method}
\alias{xsegment,TrioSetList,ANY-method}
\alias{xyplot,formula,TrioSetList-method}

\title{Class \code{"TrioSetList"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{

   A list of \code{\linkS4class{TrioSet}} objects.  Each element in the
 list is a \code{\linkS4class{TrioSet}} object for one
 chromosome. Currently, only autosomes are supported by this package, so
 instances of \code{\linkS4class{TrioSet}} typically have 22 elements.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
    \item{\code{pedigree}:}{Object of class \code{"Pedigree"} ~~ }
    \item{\code{sampleSheet}:}{Object of class \code{"SampleSheet"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
Class \code{"\linkS4class{list_or_ffdf}"}, by class "list", distance 2.
Class \code{"vectorORfactor"}, by class "list", distance 3.
}
\section{Methods}{
  \describe{
    \item{allNames}{\code{signature(object = "TrioSetList")}: ... }
    \item{[}{\code{signature(x = "TrioSetList")}: ... }
    \item{calculateMindist}{\code{signature(object = "TrioSetList")}: ... }
    \item{computeBayesFactor}{\code{signature(object = "TrioSetList")}: ... }
    \item{fatherNames}{\code{signature(object = "TrioSetList")}: ... }
    \item{initialize}{\code{signature(.Object = "TrioSetList")}: ... }
    \item{mad.marker<-}{\code{signature(x = "TrioSetList", value = "matrix")}: ... }
    \item{mad.mindist<-}{\code{signature(x = "TrioSetList")}: ... }
    \item{mad.sample<-}{\code{signature(x = "TrioSetList", value = "list")}: ... }
    \item{mindist<-}{\code{signature(object = "TrioSetList")}: ... }
    \item{mindist}{\code{signature(object = "TrioSetList")}: ... }
    \item{minimumDistance}{\code{signature(object = "TrioSetList")}: ... }
    \item{motherNames}{\code{signature(object = "TrioSetList")}: ... }
    \item{ncol}{\code{signature(x = "TrioSetList")}: ... }
    \item{nrow}{\code{signature(x = "TrioSetList")}: ... }
    \item{offspringNames}{\code{signature(object = "TrioSetList")}: ... }
    \item{pedigree}{\code{signature(object = "TrioSetList")}: ... }
    \item{prune}{\code{signature(object = "TrioSetList", ranges = "RangedDataCNV")}: ... }
    \item{sampleSheet}{\code{signature(object = "TrioSetList")}: ... }
    \item{stack}{\code{signature(object = "TrioSetList")}: creates a
 \code{TrioSet} object from a \code{TrioSetList} object}
    \item{trios}{\code{signature(object = "TrioSetList")}: ... }
    \item{xsegment}{\code{signature(object = "TrioSetList", pedigreeData = "ANY")}: ... }
    \item{xyplot}{\code{signature(x = "formula", data = "TrioSetList")}: ... }
	 }
}

\author{
R. Scharpf
}
\seealso{
  \code{\linkS4class{TrioSet}}, \code{\linkS4class{Pedigree}}

}
\examples{
showClass("TrioSetList")
data(trioSetListExample)
pedigree(trioSetList)
sampleSheet(trioSetList)
trioSet <- stack(trioSetList)
}
\keyword{classes}
