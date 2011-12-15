\name{SampleSheet}
\alias{SampleSheet}

\title{Construct an object of class SampleSheet}

\description{Construct an object of class SampleSheet}

\usage{
SampleSheet( ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%  \item{filename}{character. Ideally, the filename should be the name of
%  the file on disk.}
%  \item{id}{character.  Can be the same as filename or any unique
%  identifier for a sample. This is useful if the filenames are very long.}
%  \item{plate}{character. Chemistry plate upon which sample was processed.}
%  \item{\dots}{Additional arguments passed to \code{data.frame}}
    \item{\dots}{Additional arguments passed to \code{DataFrame}}
}

\value{An object of class \code{SampleSheet}.}

\seealso{
  \code{\linkS4class{SampleSheet}}, \code{\link{samplesheetExample}},
  \code{\linkS4class{DataFrame}}
}
\examples{
path <- system.file("extdata", package="MinimumDistance")
load(file.path(path, "sample.sheet.rda"))
samplesheetExample <- SampleSheet(filename=sample.sheet$Sample.Name,
				  id=paste("NA",substr(sample.sheet$Sample.Name, 6, 10),sep=""),
				  plate=sample.sheet$Plate,
				  well=sample.sheet$Well,
				  gender=sample.sheet$Gender,
				  SentrixPosition=sample.sheet$SentrixPosition)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}

