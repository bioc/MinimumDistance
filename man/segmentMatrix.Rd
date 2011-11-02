\name{segmentMatrix}
\alias{segment2}

\title{CBS segmentation of a minimum distance matrix or log R ratio matrix}
\description{CBS segmentation of a minimum distance matrix or log R ratio matrix}

\usage{
segment2(object, pos, chrom, verbose=TRUE, DNAcopy.verbose = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{matrix}
  \item{pos}{physical position}
  \item{chrom}{chromosome}
  \item{verbose}{Logical.  Displays progress bar if true.}
  \item{DNAcopy.verbose}{level of verbosity for DNAcopy's segment }
  \item{\dots}{additional arguments to DNAcopy's segment}
}
\details{

}
\value{\code{RangedDataCBS} object}

\seealso{
  \code{\link{xsegment}}
}
\examples{
data(trioSetListExample)
logRratios <- lrr(trioSetList[[1]])[, , 3]
##trace(segmentMatrix, browser)
rd <- segment2(logRratios, pos=position(trioSetList[[1]]), chrom=chromosome(trioSetList[[1]]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}

