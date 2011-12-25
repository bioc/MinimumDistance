\name{Pedigree}
\alias{Pedigree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Constructor for Pedigree class
}
\description{
  Constructor for Pedigree class
}
\usage{
Pedigree()
}

\details{
  Constructor for object that contains information regarding the
  case-parent trios
}
\value{
Object of class \code{Pedigree}
}

\author{
R. Scharpf
}
\seealso{
\code{\linkS4class{Pedigree}}
}
\examples{
Pedigree()
Pedigree(fatherIds=letters[1], motherIds=letters[2], offspringIds=letters[3])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classes}

