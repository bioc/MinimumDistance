\name{xsegment}
\alias{xsegment}

\title{Segment data by circular binary segmentation}
\description{
  Segment the minimum distance or the log R ratios by circular binary segmentation
}
\usage{
xsegment(object, pedigreeData, id, segment.mindist = TRUE, ..., verbose = FALSE, DNAcopy.verbose = 0)
}
\arguments{
  \item{object}{Object of class \code{TrioSetList} or \code{TrioSet}}
  \item{pedigreeData}{Object of class \code{Pedigree}}
  \item{id}{character or missing.  Name of trio to segment. If missing,
  all trios are segmented. }
  \item{segment.mindist}{Logical.  If TRUE, the minimum distance is
  segmented.  If FALSE, the log R ratios for each individual in the trio
  'id' is segmented.}
  \item{\dots}{Additional arguments are passed to the function
  \code{segment} in the DNAcopy package.}
  \item{verbose}{Logical.}
  \item{DNAcopy.verbose}{Integer. Level of verbosity for segmentation.
  See \code{segment}.}
}
\value{
  A \code{RangedDataCBS} object.
}
\author{
Rob Scharpf
}
\note{
  \code{pedigreeData} can be missing if object is of class \code{TrioSetList}.
  }

\seealso{
\code{\link[DNAcopy]{segment}}
}
\examples{
data(trioSetListExample)
##  segment the minimum distance for the first 2 chromosomes
trioSetList <- trioSetList[1:2]
minimum.distance <- calculateMindist(trioSetList)
mindist(trioSetList) <- minimum.distance
mindist.segments <- xsegment(trioSetList,
			     pedigreeData=pedigree(trioSetList),
			     id=sampleNames(trioSetList))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
