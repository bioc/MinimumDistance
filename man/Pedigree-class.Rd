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
\alias{[,Pedigree-method}
%\alias{xsegment,TrioSet,Pedigree-method}

\title{Container for storing pedigree information}
\description{
  Container for storing familial information for father, mother,
  offspring trios.
}

\section{Objects from the Class}{

  See the examples for instantiating a \code{Pedigree} object.

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

\section{Accessors}{
  \describe{

    In the following code descriptions, \code{object} is an instance of
    the \code{Pedigree} class:

    \item{}{\code{allNames(object)}: character vector of sample ids.
    Note that the length of this vector is not necessarily a multiple of
    3 if some parents have multiple offspring.}

  \item{}{\code{fatherNames(object)}: character vector of father ids (not necessarily unique).}

  \item{}{\code{motherNames(object)}: character vector of mother ids
  (not necessarily unique).}

  \item{}{\code{offspringNames(object)}: character vector of offspring
    ids (must be unique).}

  \item{}{\code{trios(object)}: \code{data.frame} of sample ids. Each
  row in the \code{data.frame} contains the ids of father, mother, and
  offspring, respectively.}

   \item{}{\code{object[i, ]}: subset the Pedigree by trio index \code{i}.}

  }

}

\author{
R. Scharpf
}

\seealso{
  \code{\link{pedigreeExample}}
}
\examples{
showClass("Pedigree")
##--------------------------------------------------
##
## instantiating an object of class Pedigree
##
##--------------------------------------------------
path <- system.file("extdata", package="MinimumDistance")
load(file.path(path, "pedigreeInfo.rda"))
trioFactor <- matrix(seq_len(nrow(pedigreeInfo)), nrow(pedigreeInfo), 3,byrow=FALSE)
indId <- unlist(pedigreeInfo)
memberId <- substr(names(indId), 1, 1)
## which row in the pedigree file
trio.index <- as.integer(trioFactor) ## note this gives the index in the pedigree
pedigreeIndex <- data.frame(individualId=indId,
			    memberId=memberId,
			    index.in.pedigree=trio.index,
			    stringsAsFactors=FALSE)
rownames(pedigreeIndex) <- NULL
pedigreeExample <- new("Pedigree", trios=pedigreeInfo,
		       trioIndex=pedigreeIndex)


##--------------------------------------------------
##
## accessors
##
##--------------------------------------------------
## Note sampleNames or offspringNames return the same thing
## as the offspring identifier uniquely identifies a trio
identical(sampleNames(pedigreeExample), offspringNames(pedigreeExample))
## To list all samples, use the accessor
allNames(pedigreeExample)
## see also
fatherNames(pedigreeExample)
motherNames(pedigreeExample)
trios(pedigreeExample)
}

\keyword{classes}
