\name{Pedigree}
\alias{Pedigree}
\title{
  Constructor for Pedigree class
}
\description{
  Constructor for Pedigree class
}
\usage{
Pedigree(pedigreeInfo, fatherIds=character(), motherIds=character(), offspringIds=character())
}

\arguments{

	\item{pedigreeInfo}{ a \code{data.frame} with column labels F, M, and O containing sample identifiers for the father (F), mother (M), and offspring (O).  Each row of the \code{data.frame} corresponds to a single trio. }

	\item{fatherIds}{ A vector of sample ids for the father.  Ignored if \code{pedigreeInfo} is specified.  }

	\item{motherIds}{ A vector of sample ids for the mother.  Ignored if \code{pedigreeInfo} is specified.  }

	\item{offspringIds}{ A vector of sample ids for the mother.  Ignored if \code{pedigreeInfo} is specified.  }

}

\details{
  Constructor for object that contains information regarding the
  case-parent trios
}
\value{
Object of class \code{Pedigree}
}

\section{Note}{

	If \code{pedigreeInfo} is missing, the vector of character
	strings provided to indicate the ids of the father, mother and
	offspring should be ordered such that the ith element in each
	vector are the sample identifiers for one trio.
}

\author{
R. Scharpf
}

\seealso{
	\code{\linkS4class{Pedigree}}
}

\examples{
	Pedigree()
	Pedigree(fatherIds=letters[1], motherIds=letters[2],offspringIds=letters[3])
	path <- system.file("extdata", package="MinimumDistance")
	load(file.path(path, "pedigreeInfo.rda"))
	Pedigree(pedigreeInfo)
}

\keyword{manip}

