\name{trioSetListExample}
\alias{trioSetListExample}
\alias{trioSetList}
\docType{data}
\title{An example of a \code{TrioSetList} object}
\description{A \code{TrioSetList} object instantiated from HapMap
  samples arrayed on a high-throughput Illumina genotyping platform.
}
\usage{data(trioSetListExample)}
\format{\code{TrioSetList} object}
\details{
  Each element in the \code{TrioSetList} is a \code{TrioSet}.
}
\source{
  Two HapMap trios.
}
\examples{
data(trioSetListExample)
class(trioSetList[[1]])
}
\keyword{datasets}
