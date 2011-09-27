\name{phenoDataArray}
\alias{phenoDataArray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an array to store the pheno data}
\description{This function creates a container to store the pheno data
for each trio based on the information provided in the pedigree and
samplesheet files.}
\usage{
phenoDataArray(pedigree, samplesheet, mapFunction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedigree}{a pedigree file}
  \item{samplesheet}{a samplesheet file}
  \item{mapFunction}{user defined function which provides a mapping of
  Sample.Name in samplesheet file to the trio identifiers in pedigree
  file. See details for more information.}
}
\details{The mapFunction is a user defined function which basically
adds a column called "pedigree.Name" in the samplesheet file which
contains a mapping between "Sample.Name" in samplesheet file and the
trio identifiers in the pedigree file.

If Sample.Name in samplesheet file and trio identifiers in
pedigree file are identical then the mapFunction is not
required. Though the user must ensure that there is a column
"pedigree.Name" in the samplesheet file which contains the same
information as in "Sample.Name".}
\value{An array with dimensions number of trios X number of columns in
samplesheet X 3. Each row represents the phenoData for one trio.}

\author{
Moiz Bootwalla
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	Pedigree files must be organized as in \code{\link{pedigreeExample}}.
	See \code{\link{samplesheetExample}} for an example samplesheet.
}
\examples{
	data(pedigreeExample)
	data(samplesheetExample)

	##Example mapFunction	
	mapFunction <- function(samplesheet) {
		pedigree.Name <- sapply(samplesheet$Sample.Name, substring, 6, 10)
		pedigree.Name <- paste("NA", pedigree.Name, sep="")
		samplesheet$pedigree.Name <- pedigree.Name
		return(samplesheet)
	}

	pData <- phenoDataArray(pedigreeExample, samplesheetExample, mapFunction)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
