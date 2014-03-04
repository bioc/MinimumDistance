\name{narrowRanges}
\alias{narrowRanges}
\title{Adjust breakpoints from segmentation inward}
\description{
Narrow the minimum distance ranges by the segmentation of the offspring
  copy number estimates.
}
\usage{
narrowRanges(object, lrr.segs, thr, mad.minimumdistance, verbose = TRUE, fD, genome)
}
\arguments{
  \item{object}{
    A \code{GRanges} object.  The segmentation of the minimum
    distance.
    }

  \item{lrr.segs}{
    A \code{GRanges} object. The segmentation of the  log R ratios.
   }

  \item{thr}{ Numeric.     For segments with a segment mean less than
    \code{thr}, the breakpoints are not altered.
  }

  \item{mad.minimumdistance}{A named numeric vector. The names should
  the \code{sampleNames} of the \code{TrioSet} or \code{TrioSetList}
  object.}

  \item{verbose}{ Logical. Whether to display messages that indicate
    progress.}

  \item{fD}{If \code{object} contains genomic ranges for a single
    chromosome, \code{fD} should be an object of class

  \code{GenomeAnnotatedDataFrame}.  If \code{object} contains genomic
  intervals from multiple chromosomes, \code{fD} should be a list of
  \code{GenomeAnnotatedDataFrame}s.  To extract a list of
  \code{GenomeAnnotatedDataFrames}, use the method \code{featureData}
  defined for a \code{TrioSetList} class. To extract a
  \code{GenomeAnnotatedDataFrame}, apply the \code{featureData} method
  to a \code{TrioSet} object.  The feature-level annotation contained in
  this object is used to update the coverage (column 'numberProbes') in the
  \code{GRanges} object after 'narrowing' the ranges.
}

\item{genome}{

  If missing, the UCSC genome build is extracted from the metadata of
  argument \code{object}.
  }

}

\details{

  If the start and stop coordinates for a segment [x, y] with mean log R
  ratio greater than \code{thr} in absolute value, the x and y
  coordinates of the interval may be adjusted.  If there are no
  breakpoints from the segmentation of the offspring log R ratio
  occurring in [x, y], nothing is done.  However, if one or more
  breakpoints occur in the interval [x,y], one or more new segments can
  be created.  For example, suppose a segment from the log R ratio
  segmention has breakpoints given by [a, b], where x < b < y.  Then the
  following two intervals are created:

  1. [x, b]
  2. [b, y]

  The motivation is to avoid having a single minimum distance segment spanning
  differing copy number states in the offspring.

}
\value{
  A  \code{GRanges} object.
}
\author{
Rob Scharpf
}

\examples{
	library(oligoClasses)
	data(trioSetListExample)
	data(lrr.segs)
	data(md.segs)
        md <- calculateMindist(lrr(trioSetList))
	md.mads <- mad2(md, byrow=FALSE)
	md.segs.narrowed <- narrowRanges(object=md.segs, lrr.segs=lrr.segs, thr=0.9, mad.minimumdistance=md.mads, fD=featureData(trioSetList))
}
\keyword{manip}

