\name{minimumDistanceCalls}
\alias{minimumDistanceCalls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Infer trio copy number state for each segment}
\description{Estimate the maximum a posteriori probability for each
  minimum distance segment.
}
\usage{
minimumDistanceCalls(id, container, chromosomes = 1:22, ranges, cbs.segs, segment.md = missing(cbs.filename) & missing(ranges), mindistance.threshold = 0.075, narrowRanges = TRUE, prOutlier = c(0.01, 1e-15), prMosaic = 0.01, prob.nonMendelian = 1.5e-06, mu.logr = c(-2, -0.5, 0, 0.3, 0.75), baf.sds = c(0.005, 0.2, 0.005), pruneMinimumDistance = FALSE, prune.by.call = TRUE, min.coverage = 10, returnEmission = FALSE, ..., verbose = TRUE, DNAcopy.verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{character.  Name of trio.}
  \item{container}{A \code{TrioSetList} object.}
  \item{chromosomes}{integer}
  \item{ranges}{A \code{RangedDataCBS} object containing the segmentation results of the minimum distance ranges.}
%  \item{cbs.filename}{filename containing CBS results}
  \item{cbs.segs}{A \code{RangedDataCBS} object.}
  \item{segment.md}{ Logical.  Whether to segment the minimum distance.}
  \item{mindistance.threshold}{Numeric. See \code{\link{narrow}}.}
  \item{narrowRanges}{Logical.  The breakpoint
    defining the start and end of a minimum distance segment are
    adjusted, if necessary, if a breakpoint from the segmentation of the
  offspring log R ratios occurs in the interval.}
  \item{prOutlier}{Numeric vector of length 2 with values in
    [0,1]. The two values are the probabilities that the log R ratio and B allele
    frequency at a marker are outliers, respectively.}
  \item{prMosaic}{Numeric.  The probability of mosaicism.}
  \item{prob.nonMendelian}{Numeric. The probability of non-Mendelian transmission.}
  \item{mu.logr}{Numeric.  The log R ratios for the copy number states.}
  \item{baf.sds}{Numeric. The standard deviation of the B allele
    frequencies for AA, AB, and BB genotypes.}
  \item{pruneMinimumDistance}{Logical. If TRUE, adjacent segments are
    concatenated (breakpoint removed) if the coverage of the 5-prime segment
    is below a threshold value given by \code{min.coverage} or if the standardized difference in
    segment means is small.}
  \item{prune.by.call}{Adjacent segments with the same trio copy number
    state call are concatenated, forming a single segment.}
  \item{min.coverage}{Numeric.}
  \item{returnEmission}{Logical.  If TRUE, an object containing the
    emission probabilities for the log R ratios and B allele frequecies
    is returned.  This is primarily for debugging -- users should set
    this value to FALSE.}
  \item{\dots}{Additional arguments passed to \code{\link{segment}}.}
  \item{verbose}{Logical.}
  \item{DNAcopy.verbose}{Integer.}
}

\details{

If \code{id} is missing, all trios are processed.

The object \code{ranges} is a \code{RangedDataCBS} object containing
the results from the segmentation of the minimum distance. If
\code{ranges} is missing, the minimum distance is segmented by the CBS
algorithm and the resulting segmentation is returned. To call the trio
copy number state for the segments, the \code{minimumDistanceCalls}
function must be called a second time passing the minimum distance
ranges to the \code{ranges} argument.  A second \code{RangedDataCBS}
object that contains the results from the segmentation of the
offspring log R ratios can be passed through the argument
\code{cbs.segs}.  If \code{narrowRanges} is TRUE, the breakpoints in
\code{ranges} for nonzero segments are adjusted such that start and
end coordinates do not span breakpoints in the \code{cbs.segs} object.

Provided that \code{ranges} is not missing, posterior probabilities
for each of the 125 possible trio copy number states are estimated for
each segment.  The inferred copy number is the state with the maximum
a posteriori probability.  If \code{prune.by.call} is true, adjacent
segments that have the same copy number state are concatenated into a
single segment.  The maximum a posteriori probability as well as the
posterior probability for the normal copy number state are returned
for each segment.

The copy number states is represented by an integer (1 + inferred copy
number).  For example, '1' corresponds to copy number zero, '2'
corresponds to copy number 1, and so on.  The copy number state for
the father, mother, and offspring in a trio is given by 3 digits.  For
example, '332' corresponds to normal copy number in the father and
mother and 1 hemizygous deletion in the offspring.  We regard copy
number alterations in the offspring that are not present in either
parent as de novo, irrespective of whether the alteration is
Mendelian.  For example, we regard the states '332', '331', '321',
'334', etc as de novo copy number states, whereas '333', '322', '311',
etc are not de novo.

}

\value{
  A \code{RangedDataCNV} object.
}

\author{
Rob Scharpf
}

\seealso{

	See \code{\link{segment}} for additional details regarding the
  segmentation provided by CBS.

  See \code{\link{gridlayout}} for plotting the low-level summaries
  for inferred deletions and amplifications.

  }

\examples{
	data(trioSetList)
	mdRanges <- minimumDistanceCalls(container=trioSetList, chromosomes=1:2)
	colnames(mdRanges)
	mdRanges <- minimumDistanceCalls(container=trioSetList,
					 chromosomes=1:2,
					 ranges=mdRanges,
					 narrowRanges=FALSE)
	colnames(mdRanges)
	state(mdRanges)
}
\keyword{models}
