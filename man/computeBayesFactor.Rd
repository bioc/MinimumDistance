\name{computeBayesFactor}
\alias{computeBayesFactor}
\alias{computeBayesFactor,TrioSet,RangedDataCNV-method}
\alias{computeBayesFactor,TrioSetList,RangedDataCNV-method}

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

  \item{\dots}{see details}
}

\details{

  Additional arguments can be specified for the non-exported functions
  \code{joint4} and \code{joint1}.  The ability to change the
  arguments for these functions from their current defaults was
  primarily useful for the development of this package. Most users
  should do fine with the default settings.

}

\value{
  A object of class \code{RangedDataCBS}.
}

\references{

  Wang et al, 2008 Nucleic Acids Research

}

\seealso{
  \code{\linkS4class{RangedDataCBS}}, \code{\link{calculateMindist}}
}

\examples{
\dontrun{
	data(trioSetListExample)
	mads.lrr.sample <- mad2(lrr(trioSetList), byrow=FALSE)
	mads.lrr.marker <- mad2(lrr(trioSetList), byrow=TRUE)
	mad.sample(trioSetList) <- mads.lrr.sample
	mad.marker(trioSetList) <- mads.lrr.marker
	mindist(trioSetList) <- calculateMindist(object=trioSetList)
	mads.md <- mad2(mindist(trioSetList), byrow=FALSE)
	md.segs <- segment2(mindist(trioSetList), chrom=chromosome(trioSetList),
			    pos=position(trioSetList))
	\dontshow{
		save(md.segs, file="~/Software/MinimumDistance/data/md.segs.rda")
	}
	lrr.segs <- segment2(lrr(trioSetList), chrom=chromosome(trioSetList),
			     pos=position(trioSetList),
			     id=trios(trioSetList))
	\dontshow{
		save(lrr.segs, file="~/Software/MinimumDistance/data/lrr.segs.rda")
	}
}
data(trioSetListExample)
data(md.segs)
data(lrr.segs)
mads.md <- mad2(mindist(trioSetList), byrow=FALSE)
md.segs2 <- narrow(md.segs, lrr.segs, mad.minimumdistance=mads.md,
		   thr=0.1)
## segmentation with maximum a posteriori probabilities (map)
map.segs <- computeBayesFactor(object=trioSetList, ranges=md.segs2)
\dontshow{
	\dontrun{
		save(map.segs, file="~/Software/MinimumDistance/man/map.segs.rda")
	}
}
}

\keyword{htest}

