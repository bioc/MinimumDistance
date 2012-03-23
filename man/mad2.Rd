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
    \code{object} can be any of the following: class \code{TrioSetList},
  \code{TrioSet}, or a list of matrices.
  }

  \item{byrow}{When \code{byrow} is \code{TRUE}, the MAD is calculated
  for each marker across all offspring. When \code{byrow} is
  \code{FALSE}, the MAD is computed across all autosomal markers for the
  fathers, mothers, and offspring. The former provides a robust estimate
  of the marker-level variance across samples, whereas the latter
  provides an estimate of variance for the samples.}

  \item{\dots}{Not currently implemented.}

}

\value{
  A list.
}

\examples{
	data(trioSetListExample)
	foreach::registerDoSEQ()
	## computing the MAD of the log R ratios for each individual (across markers)
	mads.sample <- mad2(trioSetList, byrow=FALSE)
	## compute the MAD of the log R ratio for each marker (across individuals)
	## too few samples
	##mads.marker <- mad2(trioSetList, byrow=TRUE)
}
\keyword{univar}
\keyword{robust}
