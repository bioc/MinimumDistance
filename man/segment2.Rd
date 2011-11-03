\name{segment2}
\alias{segment2}

\title{CBS segmentation of a minimum distance matrix or log R ratio matrix}
\description{CBS segmentation of a minimum distance matrix or log R ratio matrix}

\usage{
segment2(object, pos, chrom, id, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{matrix}
  \item{pos}{physical position}
  \item{chrom}{chromosome}
  \item{id}{data.frame containing trio ids.}
  \item{\dots}{additional arguments to segment function in the DNAcopy package}
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

data(trioSetListExample)
## segment the log R ratios for father mother and offspring
rd <- segment2(lrr(trioSetList),
	       pos=position(trioSetList),
	       chrom=chromosome(trioSetList),
	       id=trios(trioSetList))

## segment the minimum distance
if(is.null(mindist(trioSetList))){
	mindist(trioSetList) <- calculateMindist(trioSetList)
}
rd <- segment2(mindist(trioSetList),
	       pos=position(trioSetList),
	       chrom=chromosome(trioSetList),
	       id=trios(trioSetList))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}

