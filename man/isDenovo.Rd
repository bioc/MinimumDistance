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
  Each segment is assigned a numeric trio copy number state.  The state
  symbols are numeric.  In particular, the state "xyz" indicates that
  the father has x copies, the mother has y copies, and the offspring
  has z copies.  For example, the state '221' indicates copy number 2 in
  the mother and father, and a hemizygous deletion in the offspring.
  See \code{MAP} for details.  A subset of the possible
  state combinations, such as '221', are considered de novo.

}
\value{
Logical
}
\author{
Rob Scharpf
}

\seealso{

 \code{\link{MAP}}.  See \code{\link{state}} for the
 accessor to extract the trio copy number state.

}
\examples{
library(oligoClasses)
isDenovo(c(222, 221, 211, 223, 233))
data(map.segs)
table(isDenovo(state(map.segs)))
sts <- unique(state(map.segs))
sts[isDenovo(sts)]
}
\keyword{manip}
% Sam I am
