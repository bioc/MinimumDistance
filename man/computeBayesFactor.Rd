\name{MAP}
\alias{MAP}

\title{Estimate the trio copy number for minimum distance segments}

\description{ Estimates the trio copy number state for the minimum
  distance segments as the state that maximizes the posterior probability. The
  posterior probability for the normal state ('222') is also returned.
}

\usage{
	MAP(object, ranges, id, TAUP=1e10, tauMAX=1-5e-8, cnStates=c(-2,-0.4, 0, 0, 0.4, 1), pr.nonmendelian=1.5e-6, mdThr=0.9, ...)
}

\arguments{
  \item{object}{
    A object of class \code{TrioSetList} or \code{TrioSet}.
}

\item{ranges}{ A \code{GRanges} object. }

\item{id}{ A vector of trio identifiers (character string) for which to
  estimate the trio copy number states.  If missing, the maximum a
  posteriori probability (MAP) estimate is computed for all trios in the
  \link{TrioSetList} object.}

\item{TAUP}{ Scalar for the probability that the state of segment t is
  the same as the state at segment t-1.  Larger values decrease the
  probability of transitioning to an different state.}

\item{tauMAX}{Upper bound for the probability that the state at marker t
  is the same as the state at marker t-1.
}

\item{cnStates}{ Initial values for the mean of the copy number states
  for autosomal chromosomes.  Required to be a vector of length 6 with
  means for homozygous deletion, hemizygous deletion, diploid, diploid
  region with loss of heterozygosity, single copy gain, two-copy gain.}

\item{pr.nonmendelian}{ Scalar between 0 and 1.  The probability that
  the copy number state of the offspring was not inherited in a
  Mendelian fashion from the parental copy numbers.  }

\item{mdThr}{ Minimum distance segments with |mean| less than mdThr are
  not called.  These segments are assigned the trio-state '222'. }

\item{\dots}{Additional arguments  are passed to \code{viterbi2Wrapper}.}
}

\value{
  A  \code{GRanges} object.  The \code{value}  label "state" contains
  the trio copy number state denoted by xyz, where x is the copy number
  for the father, y is the copy number for the mother, and z is the copy
  number of the offspring.
}

\references{

  The posterior calling is an extension of the joint HMM described in
  Wang et al, 2008 Nucleic Acids Research.

}

\seealso{
	\code{\linkS4class{TrioSetList}}, \code{\linkS4class{TrioSet}}
	\code{\link{calculateMindist}}, \code{\link{viterbi2Wrapper}}
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
md.segs2 <- narrowRanges(md.segs, lrr.segs, mad.minimumdistance=mads.md, thr=0.9, fD=featureData(trioSetList))
map.segs <- MAP(trioSetList, ranges=md.segs2)
## the trio states
state(map.segs)
}
\keyword{htest}
