\name{RangedDataCNV-utils}
\alias{RangedDataCBS}
\alias{RangedDataHMM}
\alias{coverage}
\alias{state}
\alias{todf}
\alias{trioNames}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{

  Utility functions for RangedData extensions for storing ranged data on
  copy number variants.

}
\description{

  Mostly accessors for extracting data from \code{RangedDataCBS} and
  \code{RangedDataHMM} objects.
}
\usage{
RangedDataCBS(ranges = IRanges(), seg.mean = vector("numeric",  length(ranges)), ...)
RangedDataHMM(ranges=IRanges(), state=vector("integer",
  length(ranges)),...)
coverage(object)
state(object)
todf(object)
trioNames(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{RangedDataHMM} or a \code{RangedDataCNV} object.}

  \item{ranges}{
    an instance of \code{IRanges} class.
}
\item{seg.mean}{

  Numeric -- e.g., the mean copy number of a genomic interval

}
  \item{\dots}{
    Additional covariates that can be accessed by \code{$} method
}
}

\value{

  \code{coverage}, \code{state}, \code{trioNames} return a column
  (covariate) of the \code{RangedData} object.

  Coverage refers to the number of markers in the interval.  State is a
  method for \code{RangedDataHMM} objects and returns the copy number
  state, typically estimated from a HMM.

  \code{trioNames} returns the family (trio) for which the individual
  belongs.  A error message is returned if 'family' is not a column name
  in the \code{RangedData} object.


  \code{todf} coerces a \code{RangedData} object to a
  \code{\linkS4class{DataFrameCNV}} object (potentially useful for
  plotting).

}

\author{
R. Scharpf
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  See \code{\linkS4class{RangedData}} for
  additional details and methods for objects of the class.

  \code{\linkS4class{RangedDataHMM}}

}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

