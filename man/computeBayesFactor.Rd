\name{computeBayesFactor}
\alias{computeBayesFactor}
\alias{calculateMindist}
\alias{segment2}
\alias{mad2}
\alias{mad.sample<-}
\alias{mad.marker<-}
\alias{mindist<-}
\alias{narrow}

\title{Compute Bayes factors for minimum distance segments }

\description{ Estimate the trio copy number state for each minimum
  distance segment using the maximum a posterior probability. The
  posterior probability of the normal state is also returned.
}
\usage{
computeBayesFactor(object, ranges, returnEmission = FALSE, collapseRanges=TRUE, verbose = TRUE, ...)
calculateMindist(object, ..., verbose)
segment2(object, pos, chrom, ...)
mad2(object, byrow, ...)
mindist(object) <- value
mad.sample(object) <- value
mad.marker(object) <- value
narrow(md.range, cbs.segs, thr=0.9, mad.minimumdistance, verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A object of class \code{TrioSetList} or \code{TrioSet}.
}
  \item{ranges}{
    A \code{RangedDataCBS} object.
}
\item{returnEmission}{
  Whether to return the emission probabilities instead of a
  \code{RangedDataCBS} object.  This is used by the author for debugging
  and not intended for public consumption.
}

\item{collapseRanges}{ Logical.  If true, adjacent segments that have
  the same trio copy number state call are collapsed into a single segment}
\item{verbose}{Logical.}

  \item{\dots}{see details}
}

\details{

  Additional arguments can be specified for the non-exported functions
  \code{joint4} and \code{joint1}.  The ability to change the arguments
  for these functions from their current defaults was primarily useful
  for the development of this package. Most users should do fine with
  the default settings.

}

\value{
  A object of class \code{RangedDataCBS}.
}
\references{
%% ~put references to the literature/web site here ~
}

\seealso{
  \code{\linkS4class{RangedDataCBS}}, \code{\link{calculateMindist}}
}
\examples{

data(trioSetListExample)
mads.lrr.sample <- mad2(lrr(trioSetList), byrow=FALSE)
mads.lrr.marker <- mad2(lrr(trioSetList), byrow=TRUE)
mad.sample(trioSetList) <- mads.lrr.sample
mad.marker(trioSetList) <- mads.lrr.marker
mindist(trioSetList) <- calculateMindist(object=trioSetList)
mads.md <- mad2(mindist(trioSetList), byrow=FALSE)
md.segs <- segment2(mindist(trioSetList), chrom=chromosome(trioSetList),
                    pos=position(trioSetList))
lrr.segs <- segment2(lrr(trioSetList), chrom=chromosome(trioSetList),
		     pos=position(trioSetList),
		     id=trios(trioSetList))
md.segs2 <- narrow(md.segs, lrr.segs, mad.minimumdistance=mads.md,
		   thr=0.1)
rd <- computeBayesFactor(object=trioSetList, ranges=md.segs2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}

