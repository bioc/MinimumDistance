\name{xsegment}
\alias{xsegment}

\title{Segment data by circular binary segmentation}
\description{
  Segment the minimum distance or the log R ratios by circular binary segmentation
}
\usage{
xsegment(object, id, segment.mindist = TRUE, ..., verbose = FALSE, DNAcopy.verbose = 0)
}
\arguments{
  \item{object}{

}
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A \code{RangedDataCBS} object.
}
\author{
Rob Scharpf
}

\seealso{
See the function \code{\link{segment}} in the package DNAcopy.
}
\examples{
data(trioSetList)
## To segment the minimum distance for chromosome 22 -- only the first trio
mindist.segments <- xsegment(trioSetList[[22]], id=sampleNames(trioSetList)[1])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
