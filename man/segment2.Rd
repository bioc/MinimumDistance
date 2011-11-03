\name{segment2}
\alias{segment2}
\alias{segment2,array,ANY,ANY,data.frame-method}
\alias{segment2,ff_array,ANY,ANY,data.frame-method}
\alias{segment2,ff_matrix,ANY,ANY,ANY-method}
\alias{segment2,list,list,list,ANY-method}
\alias{segment2,matrix,ANY,ANY,ANY-method}

\title{Wrapper for segment function in package DNAcopy}

\description{

  \code{segment2} is a wrapper for the \code{segment} function in
  \code{DNAcopy}. The first argument can be a list, a matrix, or an
  array.  If the first argument is a list, each element of the list can
  be a matrix or an array.

  Segmentation of \code{ff_matrix} and \code{ff_array} objects is also
  supported.

  }

\usage{
segment2(object, pos, chrom, id, ...)
}

\arguments{
  \item{object}{A \code{list}, \code{matrix}, \code{array},
    \code{ff_matrix}, or \code{ff_array}.  If \code{object} is a
    \code{list}, the elements of the list must be one of the above
    classes. }

  \item{pos}{If \code{object} is a \code{list}, \code{pos} must also be
    a \code{list}.  Otherwise, \code{pos} must be an \code{integer}
    vector.}

  \item{chrom}{If \code{object} is a \code{list}, \code{chrom} must also be
    a \code{list}.  Otherwise, \code{chrom} must be an \code{integer}
    vector.}

  \item{id}{data.frame containing trio ids. The \code{id} argument is
    ignored, unless \code{object} is of class \code{array}.  See
    details.}

  \item{\dots}{Additional arguments passed to \code{segment} function in
    the \code{DNAcopy} package.}

}

\details{

  When \code{object} is a list of arrays, the argument \code{id} is
  required.  The easiest way to obtain a \code{data.frame} of the trio
  sample names is via the method \code{trios}, as in the example below.

}

\value{\code{RangedDataCBS} object}

\seealso{
	\code{\link[DNAcopy]{segment}}
      }

\examples{

data(trioSetListExample)
## segment the log R ratios for father mother and offspring
rd <- segment2(lrr(trioSetList),
	       pos=position(trioSetList),
	       chrom=chromosome(trioSetList),
	       id=trios(trioSetList))

## segment the minimum distance
if(is.null(mindist(trioSetList))){
	mindist(trioSetList) <- calculateMindist(trioSetList)
}
md.segs <- segment2(mindist(trioSetList),
	             pos=position(trioSetList),
	             chrom=chromosome(trioSetList))

}

\keyword{smooth}

