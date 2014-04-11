\name{phenoData}
\alias{fatherPhenoData}
\alias{motherPhenoData}
\title{
  Extract phenotype data for parents in a trio.
}
\description{
  Extract phenotype data for parents in a trio.
}

\usage{
motherPhenoData(object)
fatherPhenoData(object)
}

\arguments{
  \item{object}{ A object of class \code{TrioSet} or \code{TrioSetList}
  }
}
\details{
  Extracts phenotypic data for parents
}
\value{
  \code{AnnotatedDataFrame}
}

\author{
R. Scharpf
}

\seealso{
  \code{\linkS4class{TrioSet}}, \code{\linkS4class{TrioSetList}}
}

\examples{
data(trioSetListExample)
## father phenoData
fatherPhenoData(trioSetList)
## mother phenoData
motherPhenoData(trioSetList)
## offspring phenoData
phenoData(trioSetList)
}

\keyword{manip}

