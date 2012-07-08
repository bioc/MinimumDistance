\name{gcSubtract}
\alias{gcSubtract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Generic function for GC correction}
\description{Generic function for GC correction}
\usage{
gcSubtract(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A \code{TrioSet} or \code{TrioSetList} object.
}
\item{\dots}{

  Additional arguments include method and trio.index.  Method refers to
  the type of GC correction ("speed" or "lowess").  }

}

\details{

  For genomes in which the log R ratios look wavey when plotted as a
  function of physical position, a GC correction step may be helpful.
  Two options for GC correction are currently implemented. This function
  is experimental and may change in the future.

}

\value{
  An instance of the \code{object} (e.g., a \code{TrioSet} or a \code{TrioSetList}).
}

\references{
Need reference for Speed paper and reference.
}
\author{
R. Scharpf
}

\seealso{
  \code{\linkS4class{TrioSet}}, \code{\linkS4class{TrioSetList}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{smooth}
