
\name{computeBayesFactor}
\alias{computeBayesFactor}

\title{Compute Bayes' factors for minimum distance segments }

\description{ Estimate the trio copy number state for each minimum
  distance segment using the maximum a posterior probability. The
  posterior probability of the normal state is also returned.
}
\usage{
computeBayesFactor(object, ranges, returnEmission = FALSE,
  collapseRanges=TRUE, verbose = TRUE, ...)
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
\dontrun{
data(trioSetListExample)
mads.lrr.sample <- calculateMADlrr(trioSetList, by.sample=TRUE)
mads.lrr.marker <- calculateMADlrr(trioSetList, by.sample=FALSE)
mad.sample(trioSetList) <- mads.lrr.sample
mad.marker(trioSetList) <- mads.lrr.marker
md <- calculateMindist(object=trioSetList)
mads.md <- lapply(md, function(x) apply(x, 2, mad, na.rm=TRUE))
mad.mindist(trioSetList) <- mads.md
mindist(trioSetList) <- md
md.segs <- xsegment(object=trioSetList[1:10], segment.mindist=TRUE,
  verbose=TRUE)
md.segs2 <- narrow(md.segs, lrr.segs, thr=0.1)
rd <- computeBayesFactor(object=trioSetList,
			 ranges=md.segs2)

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}

