\name{RangedDataHMM-class}
\Rdversion{1.1}
\docType{class}
\alias{RangedDataHMM-class}

\title{Class \code{"RangedDataHMM"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("RangedDataHMM", ...)}.
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
Class \code{"\linkS4class{RangedDataCNV}"}, directly.
Class \code{"\linkS4class{RangedDataCopyNumber}"}, by class "RangedDataCNV", distance 2.
Class \code{"\linkS4class{RangedData}"}, by class "RangedDataCNV", distance 3.
Class \code{"\linkS4class{DataTable}"}, by class "RangedDataCNV", distance 4.
Class \code{"\linkS4class{List}"}, by class "RangedDataCNV", distance 4.
Class \code{"\linkS4class{DataTableORNULL}"}, by class "RangedDataCNV", distance 5.
Class \code{"\linkS4class{Vector}"}, by class "RangedDataCNV", distance 5.
Class \code{"\linkS4class{Annotated}"}, by class "RangedDataCNV", distance 6.
}
\section{Methods}{
No methods defined with class "RangedDataHMM" in the signature.
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
showClass("RangedDataHMM")
}
\keyword{classes}
