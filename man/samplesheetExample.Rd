\name{samplesheetExample}
\alias{samplesheetExample}
\docType{data}
\title{Container for metadata for genotyping arrays}
\description{Container for metadata for genotyping arrays}
\usage{data(samplesheetExample)}
\format{
  A data frame with 6 observations on the following 6 variables.
  \describe{
    \item{\code{filename}}{a character vector}
    \item{\code{id}}{a character vector}
    \item{\code{plate}}{a character vector}
    \item{\code{well}}{a character vector}
    \item{\code{gender}}{a character vector}
    \item{\code{SentrixPosition}}{a character vector}
  }
}

\seealso{\code{\link{SampleSheet}}, \code{\linkS4class{SampleSheet}}}

\examples{
data(samplesheetExample)
## maybe str(samplesheetExample) ; plot(samplesheetExample) ...
}
\keyword{datasets}
