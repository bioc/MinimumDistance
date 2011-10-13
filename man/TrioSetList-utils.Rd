\name{TrioSetList-utils}
\alias{allNames}
\alias{fatherNames}
\alias{fatherNames<-}
\alias{fmoNames}
\alias{motherNames}
\alias{motherNames<-}
\alias{offspringNames}
\alias{offspringNames<-}
\alias{xyplot}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{TrioSetList utility functions}

\description{ Utility functions for manipulating
  \code{\linkS4class{TrioSetList}} objects.
}

\usage{
fatherNames(object)
%fmoNames(object)
allNames(object)
motherNames(object)
offspringNames(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A \code{RangedDataCNV} object}
}
\details{
  Returns the name of the father, mother, or offspring in the trio.

%  \code{fmoNames} returns a matrix. The columns correspond to father,
%  mother, and offspring.
\code{allNames} returns all the sample names.

}
\value{
  Character string or matrix.
}
\author{
Rob Scharpf
}
\seealso{
\code{\linkS4class{TrioSet}}
}
\examples{
data("trioSetList")
fatherNames(trioSetList)
allNames(trioSetList)
motherNames(trioSetList)
offspringNames(trioSetList)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}

