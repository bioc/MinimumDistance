\name{constructTrioSetList}
\alias{constructTrioSetList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a TrioSetList object}
\description{Constructs a TrioSetList container and populates it with
logR ratios, B allele frequncies, phenoData and featureData}
\usage{
TrioSetList(pData, fD, pedigree, logR, baf, chromosome = 1:22, cdfname = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pData}{a phenoData objected created by \code{phenoDataArray}}
  \item{fD}{an object containing featureData like chromosome, position
  for each snp. See details for more information.}
  \item{pedigree}{a pedigree file}
  \item{logR}{matrix containing logR ratios for all samples}
  \item{baf}{matrix containing B allele frequencies for all samples}
  \item{chromosome}{integer}
  \item{cdfname}{character. Contains name of annotation package}
}
\details{pData must be an object returned by \code{phenoDataArray}. If
the user uses an annotation package supported by \code{Crlmm} then fD
is automatically created by the \code{TrioSetList} function. In the
above case the user must provide the correct annotation in
\code{cdfname}. Annotation packages supported by \code{Crlmm} can be
viewed using \code{annotationPackages()[grep("Crlmm",
annotationPackages())]}. If the user has annotation which is not
supported by \code{Crlmm} then they must provide the featureData.}
\value{A \code{TrioSetList} object}

\author{
Moiz Bootwalla
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	Pedigree files must be organized as in \code{\link{pedigreeExample}}.
	See \code{\link{samplesheetExample}} for an example sample sheet.
}
\examples{
\dontrun{
	data(pedigreeExample)
	data(samplesheetExample)
	data(logR)
	data(baf)
	mapFunction <- function(samplesheet) {
		pedigree.Name <- sapply(samplesheet$Sample.Name, substring, 6, 10)
		pedigree.Name <- paste("NA", pedigree.Name, sep="")
		samplesheet$pedigree.Name <- pedigree.Name
		return(samplesheet)
	}

	pData <- phenoDataArray(pedigreeExample, samplesheetExample, mapFunction)
	trioSetList <- constructTrioSetList(pData=pData, pedigree=pedigreeExample, logR=logR, baf=baf, chromosome=1:22, cdfname="human610quadv1bCrlmm")
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
