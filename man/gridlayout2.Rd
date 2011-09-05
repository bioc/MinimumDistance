\name{gridlayout2}
\alias{gridlayout2}

\title{ A visualization too for comparing the copy number breakpoints from two methods}

\description{

  This function is a wrapper for \code{grid.layout}, adding a legend
  containing posterior probabilities for the called state by
  \code{method1} and the called state and overlapping breakpoints from a
  second method.  For the second method, the breakpoints are indicated
  by blue tick marks drawn on the x-axis at the bottom.  For the first
  method, the breakpoints are indicated by orange vertical lines.

}
\usage{
gridlayout2(method1, xyList, otherCall, ranges, call.cex = 0.8, fontsize = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method1}{

    Character string.  Indicates the statistical algorithm used to call
    the ranges given in the \code{ranges} object.

  }

  \item{xyList}{The first element in the list is a list of
  \code{trellis} objects for the first range given in
  \code{ranges}. Specifically, \code{xyList[[1]][[1]]} is a
  \code{trellis} object for the copy number estimates or log R ratios
  and \code{xyList[[1]][[2]]} is a \code{trellis} object for the B
  allele frequencies.  The second element in the list corresponds to the
  second range in the \code{ranges} object and so on.




}

  \item{otherCall}{
%%     ~~Describe \code{otherCall} here~~
}
  \item{ranges}{
%%     ~~Describe \code{ranges} here~~
}
  \item{call.cex}{
%%     ~~Describe \code{call.cex} here~~
}
  \item{fontsize}{
%%     ~~Describe \code{fontsize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{aplot}

