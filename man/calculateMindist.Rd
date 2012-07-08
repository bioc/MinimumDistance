\name{calculateMindist}
\alias{calculateMindist}
\alias{calculateMindist,arrayORff_array-method}
\alias{calculateMindist,list-method}
\alias{calculateMindist,TrioSetList-method}
\title{  Compute the minimum distance. }
\description{  Compute the minimum distance.}
\usage{
calculateMindist(object, ...)
}

\arguments{

  \item{object}{ A \code{list} of arrays, a \code{TrioSetList} object,
    or an \code{array} of the log R ratios. }

  \item{\dots}{    Ignored.  }
}

\details{
  The 'minimum distance' is the minimum signed absolute difference of
    the parental log R ratios and the offspring log R
    ratios. Specifically, let |O-F| denote the absolute difference in
    the log R ratios comparing offspring to father and |O-M| the
    absolute difference in the log R ratios comparing offspring to
    mother.  The minimum distance at a marker is the signed minimum of
    |O-M| and |O-F|.  After segmentation of the minimum distance,
    non-zero segments can indicate a de novo difference in the log R
    ratio of the offspring and either parent.  For example, a positive
    minimum distance suggests that the log R ratio from the offspring is
    greater than the log R ratio of either parent.
}

\value{ If \code{object} is an array, a matrix of the minimum distance
  is returned. For an object with M markers and T trios, the dimension
  of the resulting matrix is M x T.  If \code{object} is a list of
  arrays or a \code{TrioSetList} object, a list of matrices are
  returned. The dimension of the ith matrix in the list is M_i x T.  }

\author{
R. Scharpf
}


\examples{
data(trioSetListExample)
mdlist <- calculateMindist(lrr(trioSetList))
trioSet <- stack(trioSetList)
md <- calculateMindist(lrr(trioSet))
}
\keyword{math}
\keyword{univar}
