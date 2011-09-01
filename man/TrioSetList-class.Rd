\name{TrioSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{TrioSetList-class}
\alias{[,TrioSetList-method}
\alias{fatherNames,TrioSetList-method}
\alias{motherNames,TrioSetList-method}
\alias{ncol,TrioSetList-method}
\alias{nrow,TrioSetList-method}
\alias{offspringNames,TrioSetList-method}
\alias{prune,TrioSetList,RangedDataCNV-method}
\alias{xsegment,TrioSetList-method}
\alias{xyplot,formula,TrioSetList-method}

\title{Class \code{"TrioSetList"}}
\description{


 A list of \code{\linkS4class{TrioSet}} objects.  Each element in the
 list is a \code{\linkS4class{TrioSet}} object for one
 chromosome. Currently, only autosomes are supported by this package, so
 instances of \code{\linkS4class{TrioSet}} typically have 22 elements.

}
\section{Objects from the Class}{

  Details regarding instantiation of a \code{TrioSetList} object
  are included in the vignette.

}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
Class \code{"\linkS4class{list_or_ffdf}"}, by class "list", distance 2.
Class \code{vectorORfactor}, by class "list", distance 3.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "TrioSetList")}: ... }
    \item{fatherNames}{\code{signature(object = "TrioSetList")}: ... }
    \item{motherNames}{\code{signature(object = "TrioSetList")}: ... }
    \item{ncol}{\code{signature(x = "TrioSetList")}: ... }
    \item{nrow}{\code{signature(x = "TrioSetList")}: ... }
    \item{offspringNames}{\code{signature(object = "TrioSetList")}: ... }
    \item{prune}{\code{signature(object = "TrioSetList", ranges = "RangedDataCNV")}: ... }
    \item{xsegment}{\code{signature(object = "TrioSetList")}: ... }
    \item{xyplot}{\code{signature(x = "formula", data = "TrioSetList")}: ... }
	 }
}
\author{
Rob Scharpf
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\linkS4class{TrioSet}}
}
\examples{
showClass("TrioSetList")
}
\keyword{classes}
