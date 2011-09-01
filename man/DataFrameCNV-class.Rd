\name{DataFrameCNV-class}
\Rdversion{1.1}
\docType{class}
\alias{DataFrameCNV-class}
\alias{[,DataFrameCNV-method}
\alias{xyplot,formula,DataFrameCNV-method}

\title{Class \code{"DataFrameCNV"}}
\description{
   A data.frame extension for organizing case-parent trio data.
   This class is most useful for organizing the data in a way that
   enables fairly direct application of 'xyplot' in lattice.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DataFrameCNV", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
    \item{\code{names}:}{Object of class \code{"character"} ~~ }
    \item{\code{row.names}:}{Object of class \code{"data.frameRowLabels"} ~~ }
    \item{\code{.S3Class}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{data.frame}"}, directly.
Class \code{"\linkS4class{list}"}, by class "data.frame", distance 2.
Class \code{"\linkS4class{oldClass}"}, by class "data.frame", distance 2.
Class \code{"\linkS4class{data.frameOrNULL}"}, by class "data.frame", distance 2.
Class \code{"\linkS4class{vector}"}, by class "data.frame", distance 3.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "DataFrameCNV")}: ... }
    \item{xyplot}{\code{signature(x = "formula", data = "DataFrameCNV")}: ... }
	 }
}
\author{
R. Scharpf
}
\seealso{
  \code{link{xyplot}}
}
\examples{
showClass("DataFrameCNV")
}
\keyword{classes}
