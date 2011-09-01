\name{RangedDataCNV-class}
\Rdversion{1.1}
\docType{class}
\alias{RangedDataCNV-class}
\alias{coverage,RangedDataCNV-method}
\alias{state,RangedDataCNV-method}
\alias{todf,RangedDataCNV,ANY-method}

\title{Class \code{"RangedDataCNV"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RangedDataCNV", ...)}.
%%  ~~ describe objects here ~~
}
\section{Slots}{
  \describe{
    \item{\code{ranges}:}{Object of class \code{"RangesList"} ~~ }
    \item{\code{values}:}{Object of class \code{"SplitDataFrameList"} ~~ }
    \item{\code{elementType}:}{Object of class \code{"character"} ~~ }
    \item{\code{elementMetadata}:}{Object of class \code{"DataTableORNULL"} ~~ }
    \item{\code{metadata}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
  %Class \code{"\linkS4class{RangedDataCopyNumber}"}, directly.
  Class \code{RangedDataCopyNumber}, directly.
Class \code{"\linkS4class{RangedData}"}, by class "RangedDataCopyNumber", distance 2.
Class \code{"\linkS4class{DataTable}"}, by class "RangedDataCopyNumber", distance 3.
Class \code{"\linkS4class{List}"}, by class "RangedDataCopyNumber", distance 3.
Class \code{"\linkS4class{DataTableORNULL}"}, by class "RangedDataCopyNumber", distance 4.
Class \code{"\linkS4class{Vector}"}, by class "RangedDataCopyNumber", distance 4.
Class \code{"\linkS4class{Annotated}"}, by class "RangedDataCopyNumber", distance 5.
}
\section{Methods}{
  \describe{
    \item{coverage}{\code{signature(object = "RangedDataCNV")}: ... }
    \item{prune}{\code{signature(object = "TrioSetList", ranges = "RangedDataCNV")}: ... }
    \item{prune}{\code{signature(object = "TrioSet", ranges = "RangedDataCNV")}: ... }
    \item{state}{\code{signature(object = "RangedDataCNV")}: ... }
    \item{todf}{\code{signature(object = "RangedDataCNV", range = "ANY")}: ... }
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
showClass("RangedDataCNV")
}
\keyword{classes}
