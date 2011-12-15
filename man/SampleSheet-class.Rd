\name{SampleSheet-class}
\Rdversion{1.1}
\docType{class}
\alias{SampleSheet-class}
\alias{[,SampleSheet-method}

\title{Container for storing meta-data on samples and experiment \code{"SampleSheet"}}
\description{
  Container for storing meta-data on samples and experiment
  \code{"SampleSheet"}. Future versions may use
  \code{AnnotatedDataFrame} for this.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SampleSheet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{rownames}:}{Object of class \code{"characterORNULL"} ~~ }
    \item{\code{nrows}:}{Object of class \code{"integer"} ~~ }
    \item{\code{listData}:}{Object of class \code{"list"} ~~ }
    \item{\code{elementType}:}{Object of class \code{"character"} ~~ }
    \item{\code{elementMetadata}:}{Object of class \code{"DataTableORNULL"} ~~ }
    \item{\code{metadata}:}{Object of class \code{"list"} ~~ }
  }
}

\section{Extends}{
Class \code{"\linkS4class{DataFrame}"}, directly.
Class \code{"\linkS4class{DataTable}"}, by class "DataFrame", distance 2.
Class \code{"\linkS4class{SimpleList}"}, by class "DataFrame", distance 2.
Class \code{"\linkS4class{DataTableORNULL}"}, by class "DataFrame", distance 3.
Class \code{"\linkS4class{List}"}, by class "DataFrame", distance 3.
Class \code{"\linkS4class{Vector}"}, by class "DataFrame", distance 4.
Class \code{"\linkS4class{Annotated}"}, by class "DataFrame", distance 5.
}
\section{Methods}{
  \describe{
    \item{}{\code{object[i,j,drop]}: }
  }
}

%\section{note}{ 'id' is a required column name and must be unique.}

\seealso{
  \code{\linkS4class{Pedigree}}, \code{\linkS4class{DataFrame}}
}

\examples{
showClass("SampleSheet")
## empty SampleSheet
SampleSheet()
## example
data(samplesheetExample)
}
\keyword{classes}
