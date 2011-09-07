\name{isDenovo}
\alias{isDenovo}
\title{Calculate de novo indicator}
\description{Return an indicator for each called segment for whether the
copy number of the offspring is de novo}
\usage{
isDenovo(states)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{states}{a string with 3 characters.  Character elements 1-3
    indicate the copy number state for the father, mother, and
    offspring, respectively.  See Details.
}
}
\details{
  We consider a copy number alteration occurring in the offspring that
  is not present in either parent to be de novo, irrespective of whether
  the copy number alteration is transmitted by a Mendelian mechanism.
  Each segment is assigned a trio copy number state by a 3 digit integer
  code that indicates the copy number state for the father, mother, and
  offspring.  For example, '332' indicates copy number 2 in the mother
  and father (state '3') and a hemizygous deletion in the offspring
  (state '2').    See \code{minimumDistanceCalls} for details.  A subset of the
  possible states, such as '332', is considered de novo.

}
\value{
Logical
}
\author{
Rob Scharpf
}

\seealso{
See \code{\link{minimumDistanceCalls}} for assigning the trio copy
number state to a \code{RangedDataCNV} object.  See \code{\link{state}}
for the accessor to extract the trio copy number state.
}
\examples{
isDenovo(c(333, 332, 311, 334, 344))
}
\keyword{manip}

