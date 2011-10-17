\name{calculateMADlrr}
\alias{calculateMADlrr}
%\alias{mad.marker}
\alias{mad.marker<-}
\alias{mad.sample<-}
\alias{mad.mindist<-}

\title{Compute the median absolute deviation of the log R ratio
}
\description{Compute the median absolute deviation of the log R ratio
}
\usage{
calculateMADlrr(object, by.sample)
mad.marker(x) <- value
mad.sample(x) <- value
mad.mindist(x) <- value
}

\arguments{
  \item{object}{
    Object of class \code{TrioSetList} or \code{TrioSet}
  }
 \item{x}{
    Object of class \code{TrioSetList} or \code{TrioSet}.
  }
  \item{by.sample}{Logical. If TRUE, the median absolute deviation is
  calculated across all autosomal markers for each sample.  If FALSE,
  the median absolute deviation of the log R ratios is calculated across
  samples for each marker.
}
  \item{value}{a list if object is a \code{TrioSetList}.  If object is a
  \code{TrioSet}, value must be a \code{matrix}.}
}

\details{

  \code{calculateMADlrr} calculates the median absolute deviation (MAD)
  of the log R Ratio across autosomal markers for each sample when
  \code{by.sample} is TRUE.  If \code{by.sample} is FALSE, the function
  calculates the MAD of the log R ratios for each marker across all
  offspring. The result returned by \code{calculateMADlrr} is a list of
  matrices corresponding to each element in the \code{TrioSetList}
  object.

  The replacment methods \code{mad.marker} and \code{mad.sample} are
  used to attach information on the variability of markers and samples
  to the \code{TrioSetList} object.

}

\value{
  A list.
}

\examples{
data(trioSetListExample)
## computing the MAD of the log R ratios for each individual (across markers)
mads.lrr.sample <- calculateMADlrr(trioSetList, by.sample=TRUE)
is(mads.lrr.sample, "list")
##
mad.sample(trioSetList) <- mads.lrr.sample
## compute the MAD of the log R ratio for each marker (across individuals)
mads.lrr.marker <- calculateMADlrr(trioSetList, by.sample=FALSE)
mad.marker(trioSetList) <- mads.lrr.marker
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{robust}
