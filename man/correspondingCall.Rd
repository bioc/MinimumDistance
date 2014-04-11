\name{correspondingCall}
\alias{correspondingCall}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find overlapping ranges}
\description{
  Finds ranges in a \code{RangedDataCNV} object that overlap with ranges
  from another \code{RangedDataCNV} object.
}
\usage{
correspondingCall(ranges.query, ranges.subject, subject.method)
}
\arguments{
  \item{ranges.query}{
    A \code{RangedDataCNV} object.
}
  \item{ranges.subject}{
    A \code{RangedDataCNV} object.
}
  \item{subject.method}{
    A character string.  If provided, a column called `method' will be
    added to \code{ranges.subject} that indicates the statistical
    algorithm used to call the de novo alterations.
}
}
\details{

  Given a set of de novo calls from one statistical algorithm, this
  function finds the corresponding calls made by a second statistical
  algorithm.  For any given range in \code{ranges.query}, one or more
  ranges in \code{ranges.subject} may overlap.

}
\value{
  An object of the same class as \code{ranges.subject}.
}
\author{
Rob Scharpf
}
\seealso{\code{\link{findOverlaps}}}
\keyword{utilities}

