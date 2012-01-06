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

  \item{\dots}{Additional arguments specified in the methods for \code{\linkS4class{TrioSet}} and \code{\linkS4class{TrioSetList}}.}
}

\value{
  A object of class \code{RangedDataCBS}.
}

\references{
	Wang et al, 2008 Nucleic Acids Research
}

\seealso{
	\code{\linkS4class{TrioSetList}}, \code{\linkS4class{TrioSet}}
	\code{\linkS4class{RangedDataCBS}}, \code{\link{calculateMindist}}
}

\examples{
\dontrun{
	data(trioSetListExample)
	md <- calculateMindist(lrr(trioSetList))
	md.segs <- segment2(md, chrom=as.list(chromosome(trioSetList)), pos=position(trioSetList), id=offspringNames(trioSetList), fns=featureNames(trioSetList))
	\dontshow{
		save(md.segs, file="~/Software/MinimumDistance/data/md.segs.rda")
	}
	lrr.segs <- segment2(lrr(trioSetList), chrom=as.list(chromosome(trioSetList)),
			     pos=position(trioSetList),
			     id=trios(trioSetList),
			     fns=featureNames(trioSetList))
	\dontshow{
		save(lrr.segs, file="~/Software/MinimumDistance/data/lrr.segs.rda")
	}
}
data(trioSetListExample)
data(md.segs)
data(lrr.segs)
md <- calculateMindist(lrr(trioSetList))
mads.md <- mad2(md, byrow=FALSE)
md.segs2 <- narrow(md.segs, lrr.segs, mad.minimumdistance=mads.md,
		   thr=0.1)
## segmentation with maximum a posteriori probabilities (map)
mads.lrr.sample <- mad2(lrr(trioSetList), byrow=FALSE)
mads.lrr.marker <- mad2(lrr(trioSetList), byrow=TRUE)
map.segs <- computeBayesFactor(trioSetList, ranges=md.segs2,
			       mad.marker=mads.lrr.marker,
			       mad.sample=mads.lrr.sample)
\dontshow{
	\dontrun{
		save(map.segs, file="~/Software/MinimumDistance/man/map.segs.rda")
	}
}
}

\keyword{htest}

