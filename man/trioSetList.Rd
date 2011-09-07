\name{trioSetList}
\alias{trioSetList}
\docType{data}
\title{
	Example of a TrioSetList
}
\description{

	Example of a TrioSetList. Each element of the list contains an
	object of the \code{TrioSet-class}.

}
\details{
  This object contains 22 elements, each of which is an object of the
  \code{TrioSet-class}. Each element contains the logR ratios and
  B-Allele Frequencies for each trio across one chromosome. There is
  also a slot for storing the Minimum Distance.
}
\usage{data(trioSetList)}
\examples{
data(trioSetList)
## --------------------------------------------------
## accessors for the logR Ratios and B-Allele Frequencies
## --------------------------------------------------
logR(trioSetList[[1]])[1:10, , ]
logR(trioSetList[[1]])[1:10, ,"F" ]
baf(trioSetList[[1]])[1:10, , ]
## --------------------------------------------------
## accessors for Minimum Distance calculated by a call to 
## minimumDistanceCalls
## --------------------------------------------------
mindist(trioSetList[[1]])[1:10,]
## Accessor for the MAD values calculated for each individual of a
## trio across all logR values
mad(trioSetList[[1]])
}
\keyword{datasets}
