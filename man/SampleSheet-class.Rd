\name{SampleSheet-class}
\Rdversion{1.1}
\docType{class}
\alias{SampleSheet-class}
\alias{[,SampleSheet-method}

\title{Container for storing meta-data on samples and experiment \code{"SampleSheet"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SampleSheet", ...)}.
%%  ~~ describe objects here ~~
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
    \item{[}{\code{signature(x = "SampleSheet")}: ... }
	 }
}

\seealso{
  \code{\linkS4class{Pedigree}}
}
\examples{
showClass("SampleSheet")
}
\keyword{classes}
