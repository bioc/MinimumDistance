\name{todf}
\alias{todf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coerce to data.frame}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
todf(object, range, frame, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{range}{
%%     ~~Describe \code{range} here~~
}
  \item{frame}{
%%     ~~Describe \code{frame} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (object, range, frame, ...)
standardGeneric("todf"), generic = structure("todf", package = "MinimumDistance"), package = "MinimumDistance", group = list(), valueClass = character(0), signature = c("object",
"range", "frame"), default = `NULL`, skeleton = function (object,
    range, frame, ...)
stop("invalid call in method dispatch to \"todf\" (no default method)",
    domain = NA)(object, range, frame, ...), class = structure("standardGeneric", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
