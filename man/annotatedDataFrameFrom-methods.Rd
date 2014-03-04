\name{annotatedDataFrameFrom-methods}
\docType{methods}
\alias{annotatedDataFrameFrom-methods}
\alias{annotatedDataFrameFrom,array,ANY-method}
\alias{annotatedDataFrameFrom,ff_array,ANY-method}
\alias{annotatedDataFrameFrom,ff_matrix,ANY-method}
\alias{annotatedDataFrameFrom,matrix,ANY-method}
\alias{annotatedDataFrameFrom,Pedigree,logical-method}
\title{Methods for creating AnnotatedDataFrame objects }
\description{
  Methods for creating AnnotatedDataFrame objects in the package
  \pkg{MinimumDistance}.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "array", byrow = "ANY")}}{
  Creates an AnnotatedDataFrame from an \code{array}. \code{byrow} is ignored.
}

\item{\code{signature(object = "ff_array", byrow = "ANY")}}{
  Creates an AnnotatedDataFrame from an \code{ff_array}. \code{byrow} is ignored.
}

\item{\code{signature(object = "ff_matrix", byrow = "ANY")}}{
  Creates an AnnotatedDataFrame from an \code{ff_matrix}. \code{byrow} is ignored.
}
\item{\code{signature(object = "Pedigree", byrow = "logical")}}{
    Creates an AnnotatedDataFrame from an object of class \code{Pedigree}.
}
}}

\section{Usage}{

  \describe{

    When object is of class \code{Pedigree}, an
    \code{AnnotatedDataFrame} containing phenotypic information on the
    father, mother, or offspring can be specified by the argument
    \code{which}.

    \item{}{\code{annotatedDataFrameFrom(object, byrow, sample.sheet,
	which=c("offspring", "father", "mother"), row.names=NULL, ...)}:

      A \code{data.frame} containing phenotypic information on the
      samples can be passed to this method through the argument
      \code{sample.sheet}.  The \code{sample.sheet} can contain
      phenotypic information on all the samples. Only the rows
      relevant to the offspring, for example, will be included when
      instantating an \code{AnnotatedDataFrame} when \code{which} is
      "offspring". When \code{sample.sheet} is \code{missing}, an
      \code{AnnotatedDataFrame} with zero rows will be
      instantiated. When \code{sample.sheet} is not missing,
      \code{row.names} must be specified. The \code{row.names} are the
      identifiers for each row in \code{sample.sheet} and are matched
      to sample identifiers stored in \code{object}.

    }
  }
}


\keyword{methods}

