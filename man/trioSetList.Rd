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
\details{

}
\value{
An object of class \code{TrioSetList}
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(MinimumDistance)
library(human610quadv1bCrlmm)
data(pedigreeExample)
rownames(pedigreeExample) <- NULL  ## no longer necessary to have rownames
data(samplesheetExample)
sample.sheet <- SampleSheet(filename=samplesheetExample$Sample.Name,
			    id=paste("NA",substr(samplesheetExample$Sample.Name, 6, 10),sep=""),
			    plate=samplesheetExample$Plate,
			    well=samplesheetExample$Well,
			    gender=samplesheetExample$Gender,
			    SentrixPosition=samplesheetExample$SentrixPosition)
##pedigreeTrios <- ##as(pedigreeExample, "Pedigree")
##pedigreeTrios <- as.matrix(pedigreeExample)
indId <- unlist(pedigreeExample)
trioFactor <- matrix(seq_len(nrow(pedigreeExample)), nrow(pedigreeExample), 3, byrow=FALSE)
memberId <- substr(names(indId), 1, 1)
## which row in the pedigree file
trio.index <- as.integer(trioFactor) ## note this gives the index in the pedigree
pedigreeIndex <- data.frame(individualId=indId,
			    memberId=memberId,
			    index.in.pedigree=trio.index,
			    stringsAsFactors=FALSE)
trioSetList <- TrioSetList(logR=logR,
			   baf=baf,
			   pedigree=pedigreeData,
			   sampleSheet=sample.sheet,
			   cdfname="human610quadv1bCrlmm")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}

