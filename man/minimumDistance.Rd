\name{minimumDistance}
\alias{minimumDistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the minimum distance}
\description{Calculates the minimum distance and stores the results in a
  \code{TrioSetList} container.
}
\usage{
minimumDistance(path, samplesheet, pedigree, container, container.filename, chromosomes = 1:22, cdfName, file.ext = ".txt", readFiles = TRUE, calculate.md = TRUE, calculate.mad = TRUE, exclusionRule, ..., verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{path to files containing log R ratios and B allele
    frequencies.  See details for the required format. }
  \item{samplesheet}{a \code{data.frame} providing meta-data on the samples and
    experiment.}
  \item{pedigree}{a file containing names of the trios}
  \item{container}{a \code{TrioSetList} object}
  \item{container.filename}{character string.  Filename of the
    \code{TrioSetList} object.}
  \item{chromosomes}{integer}
  \item{cdfName}{character}
  \item{file.ext}{character}
  \item{readFiles}{logical}
  \item{calculate.md}{logical.  If TRUE, the 'minimum distance' is calculated and stored in the slot 'mindist'.}
  \item{calculate.mad}{logical. If TRUE, we calculate the median absolute deviation of the log R ratios across all autosomes.  Results are stored in the slot 'mad' for each \code{TrioSet} object.}
  \item{exclusionRule}{Ignore.}
  \item{\dots}{}
  \item{verbose}{logical}
}
\details{

	Initialization of a \code{TrioSetList} object requires

	1. \code{path}: a character string indicating the complete
	path to files containing the log R ratios and B allele
	frequencies

	%% could make these requirements explicit by defining
           a formal class for the data.frame

	2. \code{samplesheet}: a \code{data.frame} object containing
	meta-data on the samples.  Required in the \code{data.frame}
	are column names 'Sample.Name', ...


	3. \code{pedigree}: a \code{matrix} indicating which samples
	belong to a trio. The sample names in the \code{pedigree}
	object should correspond to the \code{Sample.Name} column in
	the \code{samplesheet} object.

	If a \code{TrioSetList} object has already been initialized
	and is provided as an argument, \code{minimumDistance} will
	compute the signed minimum absolute difference (referred to as
	the 'minimum distance') of the log R ratios comparing the
	offspring to the parental log R ratios at each marker as well
	as the median absolute deviation (mad) of the log R ratios
	across all autosomes.  The minimum distance is stored in the
	\code{mindist} slot of each \code{TrioSet} in the
	\code{TrioSetList}.  The mad is stored in the \code{mad} slot
	of each \code{TrioSet} in the \code{TrioSetList} slot.

}
\value{A \code{TrioSetList} object.}
\author{
Rob Scharpf
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

	Pedigree files must be organized as in \code{\link{pedigreeExample}}.
	See \code{\link{samplesheetExample}} for an example sample sheet.

}
\examples{
\dontrun{
	## if a TrioSetList is provided
	data("trioSetList")
	trioSetList <- minimumDistance(container=trioSetList,
				       chromosomes=1:22,
				       cdfName="human610quadv1b",
				       calculate.md=TRUE,
				       calculate.mad=TRUE)
	save(trioSetList, file="~/Projects/MinimumDistance/data/trioSetList.rda")

	## Otherwise, the function should initialize a TrioSetList
	## object using 'path' to find files containing the BeadStudio
	## summaries, 'samplesheet' for meta-data on the samples, and
	## 'pedigree' to organize the samples by trio
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{robust}
\keyword{smooth}

