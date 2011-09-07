\name{RangedDataCBS_Examples}
\alias{cbs.segs}
\alias{md.segs}
\docType{data}
\title{Segmentation of log R ratios}
\description{
  Example \code{RangedDataCBS} objects created from the segmentation of
  the log R ratios (object \code{cbs.segs}) and the minimum distance
  (object \code{md.segs}).
}
\usage{
	data(cbs.segs)
	data(md.segs)
}
\seealso{See \code{\linkS4class{RangedDataCBS}}  for accessors and
  methods available for objects of this class. }
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
}
data(cbs.segs)
data(md.segs)
}
\keyword{datasets}
