\name{calculateMADlrr}
\alias{calculateMADlrr}

\title{Compute the median absolute deviation of the log R ratio
}
\description{Compute the median absolute deviation of the log R ratio
}
\usage{
calculateMADlrr(object, by.sample)
}

\arguments{
  \item{object}{
    Object of class \code{TrioSetList}
}
  \item{by.sample}{Logical. If TRUE, the median absolute deviation is
  calculated across all autosomal markers for each sample.  If FALSE,
  the median absolute deviation of the log R ratios is calculated across
  samples for each marker.
}
}
\details{
}
\value{
  A matrix.
}
\seealso{
  \code{\link{mad.sample}}, \code{\link{mad.marker}}
}
\examples{
data(trioSetList)
## computing the MAD of the log R ratios for each individual (across markers)
mads.lrr.sample <- calculateMADlrr(trioSetList, by.sample=TRUE)
## compute the MAD of the log R ratio for each marker (across individuals)
mads.lrr.marker <- calculateMADlrr(trioSetList, by.sample=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{robust}
