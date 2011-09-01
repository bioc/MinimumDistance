\name{pedigreeExample}
\alias{pedigreeExample}
\docType{data}
\title{
Pedigree File
}
\description{
Example of a Pedigree File. Each row contains the pedigree information
of the samples that make up a trio.
}
\usage{data(pedigreeExample)}
\format{
  A data frame with 2 observations on the following 3 variables.
  \describe{
    \item{\code{F}}{a character vector containing an identifier that
    maps back to the file containing the raw data for the Father in the trio}
    \item{\code{M}}{a character vector containing an identifier that
    maps back to the file containing the raw data for the Mother in the trio}
    \item{\code{O}}{a character vector containing an identifier that
    maps back to the file containing the raw data for the Offspring in the trio}
  }
}
\source{
HapMap Phase 2 Data
}
\examples{
data(pedigreeExample)
}
\keyword{datasets}
