\name{TrioSetList}
\alias{TrioSetList}

\title{Constructs and object of class \code{TrioSetList}}

\description{Constructs and object of class \code{TrioSetList}}
\usage{
TrioSetList(lrr, baf, pedigreeData, sampleSheet, featureData, chromosome = 1:22, cdfname)
}

\arguments{
  \item{lrr}{
    Matrix of log R Ratios
}
  \item{baf}{ matrix of B allele frequencies
}
  \item{pedigreeData}{Object of class \code{Pedigree}}
  \item{sampleSheet}{Object of class \code{SampleSheet}}
  \item{featureData}{Object of class \code{AnnotatedDataFrame}
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
library(human610quadv1bCrlmm)
data(pedigreeExample)
data(samplesheetExample)
path <- system.file("extdata", package="MinimumDistance")
load(file.path(path, "logRratio.rda"))
load(file.path(path, "baf.rda"))
trioSetList <- TrioSetList(lrr=logRratio,
			   baf=baf,
			   pedigree=pedigreeExample,
			   sampleSheet=samplesheetExample,
			   cdfname="human610quadv1bCrlmm")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}

