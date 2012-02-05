\name{mad2}
\alias{mad2}
\alias{mad2,list-method}
\alias{mad2,array-method}
\alias{mad2,matrix-method}

\title{Methods for computing the minimum absolute distance.}

\description{Compute the median absolute deviation for elements in a
	     list, array, or matrix.  For lists, elements of the list
	     can be matrices or arrays or ff-versions of these.}

\usage{
mad2(object, byrow=FALSE, ...)
}

\arguments{
  \item{object}{
    Object of class \code{TrioSetList} or \code{TrioSet}
  }

  \item{byrow}{Logical. If TRUE, the median absolute deviation is
  calculated for each row in the list elements (rows are markers).  If
  FALSE, the median absolute deviation of the log R ratios is
  calculated for each column (sample) across markers.  }

  \item{\dots}{Not currently implemented.}

}

\value{
  A list.
}

\examples{
data(trioSetListExample)
## computing the MAD of the log R ratios for each individual (across markers)
mads.sample <- mad2(lrr(trioSetList), byrow=FALSE)
## compute the MAD of the log R ratio for each marker (across individuals)
mads.marker <- mad2(lrr(trioSetList), byrow=TRUE)
}
\keyword{univar}
\keyword{robust}
