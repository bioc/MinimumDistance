\name{TrioSetList-utils}
\alias{fatherNames}
\alias{motherNames}
\alias{offspringNames}
\alias{xsegment}
\alias{xyplot}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{TrioSetList utility functions}

\description{ Utility functions for manipulating
  \code{\linkS4class{TrioSetList}} objects.
}

\usage{
fatherNames(object)
motherNames(object)
offspringNames(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{RangedDataCNV} object}
}
\details{
  Returns the name of the father, mother, or offspring in the trio.
}
\value{
  Character string
}
\author{
Rob Scharpf
}
\seealso{
\code{\linkS4class{TrioSet}}
}
\examples{
data("trioSetListExample-cleft")
fatherNames(trioSetList)
motherNames(trioSetList)
offspringNames(trioSetList)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

