\name{RangedDataCNV-utils}
\alias{RangedDataCBS}
\alias{RangedDataHMM}
\alias{RangedDataCNV}
\alias{coverage}
\alias{state}
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
RangedDataCNV(ranges = IRanges(), values, start, end, chromosome, coverage, sampleId, startIndexInChromosome, endIndexInChromosome, ...)
RangedDataCBS(ranges = IRanges(), seg.mean = vector("numeric",  length(ranges)), ...)
RangedDataHMM(ranges=IRanges(), state=vector("integer",length(ranges)),...)
coverage(object)
state(object)
trioNames(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{RangedDataHMM} or a \code{RangedDataCNV}
    object.}

  \item{ranges}{
    an instance of \code{IRanges} class.}

    \item{values}{A \code{DataFrame} object.}
    \item{start}{integer -- physical position indicating start of copy number
      segment}
    \item{end}{integer -- physical position indicating end of copy
      number segment}
    \item{chromosome}{integer indicating chromosome}
    \item{sampleId}{character string}
    \item{startIndexInChromosome}{index of marker within the chromosome
    for the beginning of the copy number segment.  The physical position
    of this marker is given by 'start'. Optional}
    \item{endIndexInChromosome}{index of marker within the chromosome
    for the end of the copy number segment.  The physical position
    of this marker is given by 'end'. Optional}
  \item{seg.mean}{
    Numeric -- e.g., the mean copy number of a genomic interval
  }
  \item{coverage}{number of markers in a segment}
  \item{state}{typically an integer corresponding to the inferred copy number from a hidden
    markov model}
  \item{\dots}{
    Additional covariates that can be accessed by \code{$} method}
}

\details{ \code{RangedDataCNV}, \code{RangedDataHMM}, and
  \code{RangedDataCBS} are constructors for the corresponding class.}

\value{

  \code{coverage}, \code{state}, \code{trioNames} return a column
  (covariate) of the \code{RangedData} object.

  Coverage refers to the number of markers in the interval.  State is a
  method for \code{RangedDataHMM} objects and returns the copy number
  state, typically estimated from a HMM.

  \code{trioNames} returns the family (trio) for which the individual
  belongs.  A error message is returned if 'family' is not a column name
  in the \code{RangedData} object.


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
  ranges <- IRanges(c(1,2,3),c(4,5,6))
  chrom <- 1:3
  id <- letters[1:3]
  num.mark <- rpois(3, 10)
  seg.mean <- rnorm(3)
  rd <- RangedDataCBS(ranges=ranges,
		      chromosome=chrom,
		      sampleId=id,
		      coverage=num.mark,
		      seg.mean=seg.mean)
  ## accessors
  chromosome(rd)
  sampleNames(rd)
  coverage(rd)

  rd.hmm <- RangedDataHMM(ranges=ranges,
			  chromosome=chrom,
			  sampleId=id,
			  num.mark=num.mark,
			  state=2L)
  ## accessors
  chromosome(rd.hmm)
  sampleNames(rd.hmm)
  coverage(rd.hmm)
  state(rd.hmm)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

