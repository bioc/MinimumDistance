\name{prune}
\alias{prune}
\alias{prune,TrioSetList,RangedDataCNV-method}

\title{Removing breakpoints from a segmented genome}

\description{

  A function for removing splits as a function of the standardized
  difference in adjacent segment means and coverage.  We use an
  exponential rule of the form
  min.change+lambda*exp(-lambda*coverage)/scale.exp. The rule is applied
  recursively until no segments are removed.

}
\usage{
prune(object, ranges, id, lambda = 0.05, min.change = 0.1, min.coverage = 3, scale.exp = 0.02, verbose, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \code{\linkS4class{TrioSet}} or \code{\linkS4class{TrioSetList}}.
}
  \item{ranges}{

    A object derived from class \code{\linkS4class{RangedDataCNV}}.

}
  \item{id}{

    offspring id

  }

  \item{lambda}{

    See the exponential curve defined in the description.

  }

  \item{min.change}{

    the minimum standardized difference in segment means

  }

  \item{min.coverage}{

    the minimum coverage

  }

  \item{scale.exp}{     See the exponential curve defined in the description.}
  \item{verbose}{ Logical -- verbose output during processing}
  \item{\dots}{Not implemented}
}

\value{
  A \code{\linkS4class{RangedDataCBS}} object.
}
\author{
Rob Scharpf
}

\seealso{
  \code{\link{segment}}
}
\keyword{manip}

