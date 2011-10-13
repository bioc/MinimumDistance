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
\alias{xsegment,TrioSet,Pedigree-method}

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
Class \code{"\linkS4class{vectorORfactor}"}, by class "list", distance 3.
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
rownames(pedigreeExample) <- NULL  ## no longer necessary to have rownames
data(samplesheetExample)
sample.sheet <- SampleSheet(filename=samplesheetExample$Sample.Name,
			    id=paste("NA",substr(samplesheetExample$Sample.Name, 6, 10),sep=""),
			    plate=samplesheetExample$Plate,
			    well=samplesheetExample$Well,
			    gender=samplesheetExample$Gender,
			    SentrixPosition=samplesheetExample$SentrixPosition)
indId <- unlist(pedigreeExample)
trioFactor <- matrix(seq_len(nrow(pedigreeExample)), nrow(pedigreeExample), 3, byrow=FALSE)
memberId <- substr(names(indId), 1, 1)
## which row in the pedigree file
trio.index <- as.integer(trioFactor) ## note this gives the index in the pedigree
pedigreeIndex <- data.frame(individualId=indId,
			    memberId=memberId,
			    index.in.pedigree=trio.index,
			    stringsAsFactors=FALSE)
rownames(pedigreeIndex) <- NULL
pedigreeData <- new("Pedigree", trios=pedigreeExample,
		    trioIndex=pedigreeIndex)
allNames(pedigreeData)
fatherNames(pedigreeData)
motherNames(pedigreeData)
offspringNames(pedigreeData)
}

\keyword{classes}
