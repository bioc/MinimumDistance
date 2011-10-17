\name{Pedigree-class}
\Rdversion{1.1}
\docType{class}
\alias{Pedigree-class}
\alias{allNames,Pedigree-method}
\alias{fatherNames,Pedigree-method}
\alias{initialize,Pedigree-method}
\alias{motherNames,Pedigree-method}
\alias{offspringNames,Pedigree-method}
\alias{trios,Pedigree-method}
%\alias{xsegment,TrioSet,Pedigree-method}

\title{Container for storing pedigree information}
\description{
  Container for storing familial information for father, mother,
  offspring trios.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Pedigree", trios, trioIndex)}.

}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} ~~ }
    \item{\code{trios}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{trioIndex}:}{Object of class \code{"data.frame"} ~~ }
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
    \item{allNames}{\code{signature(object = "Pedigree")}: ... }
    \item{fatherNames}{\code{signature(object = "Pedigree")}: ... }
    \item{initialize}{\code{signature(.Object = "Pedigree")}: ... }
    \item{motherNames}{\code{signature(object = "Pedigree")}: ... }
    \item{offspringNames}{\code{signature(object = "Pedigree")}: ... }
    \item{trios}{\code{signature(object = "Pedigree")}: ... }
    \item{xsegment}{\code{signature(object = "TrioSet", pedigreeData = "Pedigree")}: ... }
	 }
}

\author{
R. Scharpf
}

\seealso{
  \code{\linkS4class{SampleSheet}}
}
\examples{
showClass("Pedigree")
data(pedigreeExample)
## Note that
sampleNames(pedigreeExample)
## is the same as
offspringNames(pedigreeExample)
## To list all samples, use the accessor
allNames(pedigreeExample)
## see also
fatherNames(pedigreeExample)
motherNames(pedigreeExample)
}

\keyword{classes}
