\name{isDenovo}
\alias{isDenovo}
\title{Indicator for whether a trio copy number state is a de novo copy
  number alteration. }
\description{Return an indicator for a de novo copy number alteration.}
\usage{
isDenovo(states)
}
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
  (state '2').  See \code{computeBayesFactor} for details.  A subset of
  the possible states, such as '332', is considered de novo.

}
\value{
Logical
}
\author{
Rob Scharpf
}

\seealso{

 \code{\link{computeBayesFactor}}.  See \code{\link{state}} for the
 accessor to extract the trio copy number state.

}
\examples{
library(oligoClasses)
isDenovo(c(333, 332, 311, 334, 344))
data(map.segs)
table(isDenovo(state(map.segs)))
sts <- unique(state(map.segs))
sts[isDenovo(sts)]
}
\keyword{manip}
% Sam I am
