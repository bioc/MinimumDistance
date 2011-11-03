\name{mad2}
\alias{mad2}
%\alias{mad.marker}
%\alias{mad.marker<-}
%\alias{mad.sample<-}
%\alias{mad.mindist<-}

\title{Compute the median absolute deviation of the log R ratio
}
\description{Compute the median absolute deviation of the log R ratio
}
\usage{
mad2(object, by.sample)
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
mads.lrr.sample <- mad2(lrr(trioSetList), byrow=FALSE)
##
mad.sample(trioSetList) <- mads.lrr.sample
## compute the MAD of the log R ratio for each marker (across individuals)
mads.lrr.marker <- mad2(lrr(trioSetList), byrow=TRUE)
mad.marker(trioSetList) <- mads.lrr.marker
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{robust}
