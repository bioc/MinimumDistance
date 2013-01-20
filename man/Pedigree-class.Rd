\name{Pedigree-class}
\Rdversion{1.1}
\docType{class}
\alias{Pedigree-class}
\alias{allNames,Pedigree-method}
\alias{dim,Pedigree-method}
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
    \item{\code{trios}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{trioIndex}:}{Object of class \code{"data.frame"} ~~ }
  }
}

\section{Accessors}{
  \describe{

    In the following code descriptions, \code{object} is an instance of
    the \code{Pedigree} class:

    \item{}{\code{allNames(object)}: character vector of unique sample ids.
    Note that the length of this vector is not necessarily a multiple of
    3 if some parents have multiple offspring.}

  \item{}{\code{annotatedDataFrameFrom(object, byrow=FALSE, sample.sheet,
    which=c("offspring", "father", "mother"),
    row.names=NULL, ...)}:  Create an AnnotatedDataFrame for the father,
    mother or, offspring (depending on value of argument \code{which}).
    Optionally, a \code{data.frame} of covariates for the samples can be
    passed as an argument to \code{sample.sheet}.  If
    \code{sample.sheet} is not missing, \code{row.names} can not be
    \code{NULL}.}

  \item{}{\code{dim(object)}: Returns an integer vector of length 2: the
  first element is the number of trios; the second element is 3.}

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
  \code{\link{Pedigree}} \code{\link{annotatedDataFrameFrom}}
}

\examples{
showClass("Pedigree")
## an empty container
Pedigree()
}

\keyword{classes}
