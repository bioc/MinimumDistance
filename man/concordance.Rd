\name{concordance}
\alias{concAtTop}
\alias{discAtTop}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions for assessing concordance}
\description{Functions for assessing concordance and discordance of copy number variant calls}
\usage{
concAtTop(ranges.query, ranges.subject, list.size, verbose = TRUE, ...)
discAtTop(ranges.query, ranges.subject, verbose=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ranges.query}{
    \code{RangedDataCNV} object
}
  \item{ranges.subject}{
    \code{RangedDataCNV} object
}
  \item{list.size}{size of list in query and subject ranges}
  \item{verbose}{logical}
  \item{\dots}{Additional arguments passed to \code{findOverlaps}.}
}
\details{
  \code{concAtTop} calculates three measures of concordance:

  1.  the proportion of top ranges that overlap between
  \code{ranges.query} and \code{ranges.subject} objects as a function of
  list size (we assume that each \code{RangedDataCNV} object is ordered
  such that the first range has the highest rank (most evidence of an
  alteration).

  2.  the proportion of top ranges in \code{ranges.query} that appear
  anywhere in the \code{ranges.subject} object.  Again, a proportion is
  calculated as a function of list size until the size of the list is
  equal to \code{list.size}.

  3. the same as (2), but reversing the role of \code{ranges.query} and
  \code{ranges.subject}.

  The function \code{discAtTop} identifies the ranges in
  \code{ranges.query} that do not appear in \code{ranges.subject}.

}
\value{
  \code{distAtTop} returns an object of the same class as
  \code{ranges.query}.

  \code{concAtTop} returns a list of 3 elements, corresponding to the 3
  approaches for estimating concordance described in the details.

}
\author{
Rob Scharpf
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{findOverlaps}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
\keyword{manip}
