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
  \code{\link{minimumDistanceCalls}} for details regarding how the
  object \code{called.ranges} was instantiated.}

\examples{
\dontrun{
	data("trioSetList")
	## how the 'md.segs' object was created
	md.segs <- xsegment(object=trioSetList,
			    verbose=TRUE)
	## how the 'cbs.segs' object was created
	cbs.segs <- xsegment(object=trioSetList,
			     segment.mindist=FALSE,
			     verbose=TRUE)

        ## how the 'called.ranges' object was created
	called.ranges <- minimumDistanceCalls(container=trioSetList,
					      chromosomes=1:2,
					      cbs.segs=cbs.segs,
					      ranges=md.segs)


}
data(cbs.segs)
data(md.segs)
}
\keyword{datasets}
