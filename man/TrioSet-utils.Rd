\name{TrioSet-utils}
\alias{baf}
\alias{baf<-}
\alias{logR}
\alias{logR<-}
\alias{mad}
\alias{mad<-}
\alias{mindist}
\alias{mindist<-}

\title{Utility functions for TrioSet}
\description{ Utility functions for \code{TrioSet} objects}
\usage{
baf(object)
baf(object) <- value
logR(object)
logR(object) <- value
mad(x, center=median(x), constant=1.4826, na.rm=FALSE,low=FALSE,high=FALSE)
mad(object) <- value
mindist(object)
mindist(object) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{\linkS4class{TrioSet}}}
  \item{x}{Object of class \code{\linkS4class{TrioSet}}}
  \item{center}{ignored}
  \item{constant}{ignored}
  \item{na.rm}{ignored}
  \item{low}{ignored}
  \item{high}{ignored}
  \item{value}{

    An array or matrix.

  }
}

\value{
  An array or matrix
}

\author{
Rob Scharpf
}

\seealso{
\code{\linkS4class{TrioSet}}
}
\examples{
data("trioSetListExample-cleft")
trioSet <- trioSetList[[1]]
## log R ratios for the first trio and markers 1-5
logR(trioSet)[1:5, 1, ]
## log R ratios for the father for the first 5 trios and markers 1-5
logR(trioSet)[1:5, 1:5, "F"]

## the minimum distance for the first five markers, first trio
mindist(trioSet)[1:5, 1]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{manip}
