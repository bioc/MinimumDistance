\name{RangedDataCBS_Examples}
\alias{md.segs}
\alias{lrr.segs}
\alias{map.segs}
\docType{data}
\title{Objects containing segmented data}
\description{

  Example \code{RangedDataCBS} objects created from the segmentation of
  the log R ratios (object \code{cbs.segs}) and the minimum distance
  (object \code{md.segs}).

}
\usage{
	data(lrr.segs)
	data(md.segs)
        data(map.segs)
}

\seealso{

  See \code{\linkS4class{RangedDataCBS}} for accessors and methods
  available for objects of this class.  See
  \code{\link{MAP}} for details regarding how the object
  \code{map.segs} was instantiated.

}

\examples{
data(lrr.segs)
data(md.segs)


###object containing maximum a posterior probabilities
data(map.segs)
}
\keyword{datasets}
