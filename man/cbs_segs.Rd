\name{RangedDataCBS_Examples}
\alias{cbs.segs}
\alias{md.segs}
\alias{called.ranges}
\docType{data}
\title{Segmentation of log R ratios}
\description{

  Example \code{RangedDataCBS} objects created from the segmentation of
  the log R ratios (object \code{cbs.segs}) and the minimum distance
  (object \code{md.segs}).

  The object \code{called.ranges} is returned by the
  \code{minimumDistanceCalls} function and differs from \code{md.segs}
  in that a copy number state has been assigned to each segment and
  adjacent segments that have the same copy number state are combined
  into a single segment.


}
\usage{
	data(cbs.segs)
	data(md.segs)
        data(called.ranges)
}

\seealso{See \code{\linkS4class{RangedDataCBS}} for accessors and
  methods available for objects of this class.  See
  \code{\link{computeBayesFactor}} for details regarding how the
  object \code{called.ranges} was instantiated.}

\examples{
	data(called.ranges)
}
\keyword{datasets}
