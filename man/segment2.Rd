\name{segment2}
\alias{segment2}
\alias{segment2,TrioSet-method}
\alias{segment2,TrioSetList-method}
\alias{segment2,ff_matrix-method}
\alias{segment2,arrayORff_array-method}
\alias{segment2,list-method}
\alias{segment2,matrix-method}

\title{Wrapper for segment function in package DNAcopy}

\description{

  \code{segment2} is a wrapper for the \code{segment} function in
  \code{DNAcopy}. The first argument can be simple data structures: a
  \code{list}, a \code{matrix}, or an \code{array}, or more complex:
  \code{TrioSetList} and \code{TrioSet}.  If the first argument is a
  list, each element of the list can be a matrix or an array.

  Segmentation of \code{ff_matrix} and \code{ff_array} objects is also
  supported.

  }

\section{Arguments}{

    Arguments to \code{segment2} depend on the class of
    \code{object}. In all forms, additional arguments to the
    \code{segment} function in the \pkg{DNAcopy} package can be passed
    through the \code{...} operator.

    For \code{TrioSetList} objects, the arguments are:

    \describe{
      \item{}{\code{segment2(object, md=NULL, segmentParents=TRUE,
      verbose=TRUE, ...)}:   \code{md} is a list of matrices containing
    the minimum distance. See \code{\link{calculateMindist}}.
    \code{segmentParents} must be logical. When TRUE, the parental log
    R ratios are segmented via circular binary segmentation.  When
    FALSE, only the offspring log R ratios are segmented.  }

  For objects of class \code{TrioSet}:

  \item{}{\code{segment2(object, md=NULL, segmentParents=TRUE,
  verbose=TRUE, ...)}: \code{md} is a matrix of the minimum distance
  where each column corresponds to a trio in the \code{TrioSetList}
  object.}

For objects of class \code{list}:

\item{}{\code{segment2(object, pos, chrom, id=NULL, featureNames,
    segmentParents=TRUE, verbose=TRUE, ...)}: \code{pos} is a list of
    the genomic positions (integers) for each row of the elements of the
    \code{object} list.  Similarly, \code{chrom} and \code{featureNames}
    are lists specifying the chromosome (\code{integer}) and feature
    identifiers (\code{character}) for each row in the \code{object}
    list elements, respectively. Note: \code{pos}, \code{chrom}, and
    \code{featureNames} must be lists of the same size.}

  For objects of class \code{matrix} :

  \item{}{\code{segment2(object, pos, chrom, id=NULL, featureNames,
segmentParents=TRUE, verbose=TRUE, ...)}: \code{pos} is a vector of the
genomic positions (integers) for each row in the \code{object} matrix.
Similarly, \code{chrom} and \code{featureNames} are vectors specifying
the chromosome (\code{integer}) and feature identifiers
(\code{character}) for each row in the \code{object} matrix,
respectively. Note: \code{pos}, \code{chrom}, and \code{featureNames}
must be vectors of the same size.}

}
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

\dontrun{ ## examples are checked in the vignette
data(trioSetListExample)
mdlist <- calculateMindist(lrr(trioSetList))
md.segs <- segment2(trioSetList, md=mdlist)
lrr.segs <- segment2(trioSetList, segmentParents=FALSE)
}

}

\keyword{smooth}

