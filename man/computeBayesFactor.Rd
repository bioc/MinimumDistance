\name{computeBayesFactor}
\alias{computeBayesFactor}

\title{Compute Bayes factors for minimum distance segments }

\description{ Estimate the trio copy number state for each minimum
  distance segment using the maximum a posterior probability. The
  posterior probability of the normal state is also returned.
}

\usage{
	computeBayesFactor(object, ranges, ...)
}

\arguments{
  \item{object}{
    A object of class \code{TrioSetList} or \code{TrioSet}.
}

  \item{ranges}{
    A \code{RangedDataCBS} object.
}

  \item{\dots}{Additional arguments specified in the methods for
    \code{\linkS4class{TrioSet}} and \code{\linkS4class{TrioSetList}}.}

}

\value{
  A object of class \code{RangedDataCBS}.
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
data(trioSetListExample)
data(md.segs)
data(lrr.segs)
mdlist <- calculateMindist(lrr(trioSetList))
mads.md <- mad2(mdlist, byrow=FALSE)
md.segs2 <- narrow(md.segs, lrr.segs, mad.minimumdistance=mads.md, thr=0.8)
map.segs <- computeBayesFactor(trioSetList, ranges=md.segs2)
}
\keyword{htest}
