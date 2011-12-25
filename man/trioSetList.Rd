\name{TrioSetList}
\alias{TrioSetList}

\title{Constructs and object of class \code{TrioSetList}}

\description{Constructs and object of class \code{TrioSetList}}
\usage{
TrioSetList(lrr, baf, pedigreeData, sampleSheet, featureAnnotation, chromosome = 1:22, cdfname)
}

\arguments{
  \item{lrr}{
    Matrix of log R Ratios
}
  \item{baf}{ matrix of B allele frequencies
}
  \item{pedigreeData}{Object of class \code{Pedigree}}
  \item{sampleSheet}{Object of class \code{SampleSheet}}
  \item{featureAnnotation}{Object of class \code{AnnotatedDataFrame}
  containing feature annotation }
  \item{chromosome}{integer vector. Which chromosomes to include in the
  TrioSetList object.  }
\item{cdfname}{character.  Name of annotation package.}
}

\value{
An object of class \code{TrioSetList}
}

\seealso{
  \code{\linkS4class{TrioSetList}}
}
\examples{
TrioSetList()
}
\keyword{classes}

