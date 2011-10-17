
\name{trios}
\alias{trios}
\title{ Accessor for trio sample names}

\description{Accessor for the \code{trios} slot in objects of class \code{Pedigree}}
\usage{
trios(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    Object of class \code{Pedigree}, or a \code{TrioSetList}.
  }
}

\details{

  Accessor for slot \code{trios} for objects of class \code{Pedigree}

}

\value{

  Returns the sample membership in father (F), mother (M), offspring (O)
  trios (each row of the \code{data.frame} is a complete F-M-O trio).

  Parents with multiple offspring can be listed in separate rows -- the
  parental names can be duplicated, but the offspring identifier can
  not.

}

\seealso{
\code{\linkS4class{Pedigree}}
}
\examples{
data(pedigreeExample)
trios(pedigreeExample)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

