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
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TrioSetList", ...)}.
%%  ~~ describe objects here ~~
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
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("TrioSetList")
}
\keyword{classes}
