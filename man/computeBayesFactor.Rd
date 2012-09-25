\name{computeBayesFactor}
\alias{computeBayesFactor}

\title{Compute Bayes factors for minimum distance segments }

\description{ Estimates the trio copy number state for the minimum
  distance segments as the state that maximizes the posterior probability. The
  posterior probability for the normal state ('333') is also returned.
}

\usage{
	computeBayesFactor(object, ranges, ...)
}

\arguments{
  \item{object}{
    A object of class \code{TrioSetList} or \code{TrioSet}.
}

  \item{ranges}{ A \code{GRanges} object. }

  \item{\dots}{Additional arguments specified in the methods for
    \code{\linkS4class{TrioSet}} and \code{\linkS4class{TrioSetList}}.}

}

\value{
  A  \code{GRanges} object.
}

\references{

  The posterior calling is an extension of the joint HMM described in
  Wang et al, 2008 Nucleic Acids Research.

}

\seealso{
	\code{\linkS4class{TrioSetList}}, \code{\linkS4class{TrioSet}}
	\code{\linkS4class{RangedDataCBS}}, \code{\link{calculateMindist}}
}

\examples{
library(MinimumDistance)
library(oligoClasses)
foreach::registerDoSEQ()
data(trioSetListExample)
data(md.segs)
data(lrr.segs)
mdlist <- calculateMindist(lrr(trioSetList))
mads.md <- mad2(mdlist, byrow=FALSE)
md.segs2 <- narrow(md.segs, lrr.segs, mad.minimumdistance=mads.md, thr=0.9, fD=featureData(trioSetList))
map.segs <- computeBayesFactor(trioSetList, ranges=md.segs2, prOutlierBAF=list(initial=1e-4, max=1e-2, maxROH=1e-4))
}
\keyword{htest}
