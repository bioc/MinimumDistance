\name{TrioSet-utils}
%\alias{baf}
\alias{baf<-}
%\alias{lrr}
\alias{lrr<-}
\alias{mindist}
\alias{mindist<-}
\alias{phenoData2}


\title{Utility functions for TrioSet}
\description{ Utility functions for \code{TrioSet} objects}
\usage{
%baf(object)
baf(object) <- value
%lrr(object)
lrr(object) <- value
mindist(object)
mindist(object) <- value
phenoData2(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{\linkS4class{TrioSet}}}
  \item{value}{ For \code{lrr} and \code{baf} methods, \code{value} is a
    3-dimensional array.  For \code{mindist}, \code{value} is a matrix.}
}

\details{

  \code{baf} is an accessor for the B allele frequencies.

  \code{lrr} is an accessor for the log R ratios.

  \code{mindist} is an accessor for the minimum distance.

  \code{phenoData2} is an accessor for phenotypic information and
   metadata on the trios stored as an array (trios x covariate x trio
   member).

}

\value{
  An array
}

\author{
Rob Scharpf
}

\seealso{
\code{\linkS4class{TrioSet}}, \code{\link{calculateMindist}}
}
\examples{
data(trioSetListExample)
trioSet <- trioSetList[[1]]
## log R ratios for the first trio and markers 1-5
lrr(trioSet)[1:5, 1, ]
## the minimum distance for the first five markers, first trio
mindist(trioSet)[1:5, 1]
}

\keyword{manip}
