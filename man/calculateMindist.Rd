\name{calculateMindist}
\alias{calculateMindist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the minimum distance }
\description{Calculates the minimum distance }
\usage{
calculateMindist(object, ..., verbose=TRUE)
}
\arguments{
  \item{object}{
    \code{TrioSetList} or \code{TrioSet} object.}
  \item{\dots}{Ignored}
  \item{verbose}{Logical. Displays a bar indicating progress if TRUE}
}

\details{

  Let |O-F| denote the absolute difference in the log R ratios comparing
  offspring to father and |O-M| the absolute difference in the log R
  ratios comparing offspring to mother.

  The minimum distance at a marker is the signed minimum of |O-M| and
  |O-F|.  Plotted for all markers, deflections from zero can indicate a
  difference in the log R ratio of the offspring and either parent
  (possibly de novo).  For example, a positive minimum distance suggests
  that the log R ratio from the offspring is greater than the log R
  ratio of either parent.

}

\value{If \code{object} is a \code{TrioSetList}, the returned value is a
  list where each element is a matrix of the minimum distance for the
  corresponding list element in \code{TrioSetList}.

  If \code{object} is a \code{TrioSet}, the returned value is a matrix.

}

\seealso{
}
\examples{
data(trioSetListExample)
mdlist <- calculateMindist(trioSetList)
}
\keyword{univar}

