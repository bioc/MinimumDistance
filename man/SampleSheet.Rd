\name{SampleSheet}
\alias{SampleSheet}

\title{Constructs and object of class SampleSheet }

\description{Constructs and object of class SampleSheet }
\usage{
SampleSheet(filename, id, plate, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{character}
  \item{id}{character.  Can be the same as filename or any unique
  identifier for a sample.}
  \item{plate}{character. Chemistry plate upon which sample was processed.}
  \item{\dots}{Additional arguments passed to \code{data.frame}}
}
\details{

}
\value{An object of class \code{SampleSheet}.
}

\seealso{
  \code{\linkS4class{SampleSheet}}
}
\examples{
data(samplesheetExample)
sample.sheet <- SampleSheet(filename=samplesheetExample$Sample.Name,
			    id=paste("NA",substr(samplesheetExample$Sample.Name, 6, 10),sep=""),
			    plate=samplesheetExample$Plate,
			    well=samplesheetExample$Well,
			    gender=samplesheetExample$Gender,
			    SentrixPosition=samplesheetExample$SentrixPosition)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}

