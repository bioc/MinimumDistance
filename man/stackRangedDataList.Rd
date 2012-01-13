\name{stackRangedDataList}
\alias{stackRangedDataList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Stack a list of RangedDataCBS objects.
}
\description{
  Stack a list of RangedDataCBS objects.
}
\usage{
stackRangedDataList(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
    A list. Each element has class \code{RangedDataCBS}.
}
}

\value{
Object of class \code{RangedDataCBS}.
}

\author{
R. Scharpf
}

\seealso{
\code{\linkS4class{RangedDataCBS}}
}
\examples{
data(map.segs)
stackRangedDataList(list(map.segs, map.segs))
}
\keyword{manip}

