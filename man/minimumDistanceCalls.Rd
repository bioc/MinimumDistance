\name{minimumDistanceCalls}
\alias{minimumDistanceCalls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate trio copy number state for each segment}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
minimumDistanceCalls(id, container, chromosomes = 1:22, ranges, cbs.filename, cbs.segs.offspring, segment.md = missing(cbs.filename) & missing(ranges), offspring.segs, mindistance.threshold = 0.075, narrowRanges = TRUE, calculate.lr = TRUE, prOutlier = c(0.01, 1e-15), prMosaic = 0.01, prob.nonMendelian = 1.5e-06, mu.logr = c(-2, -0.5, 0, 0.3, 0.75), baf.sds = c(0.005, 0.2, 0.005), pruneMinimumDistance = FALSE, prune.by.call = TRUE, min.coverage = 10, returnEmission = FALSE, ..., verbose = TRUE, DNAcopy.verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{character.  Name of trio.}
  \item{container}{A \code{TrioSetList} object.}
  \item{chromosomes}{integer}
  \item{ranges}{}
  \item{cbs.filename}{filename containing CBS results}
  \item{cbs.segs.offspring}{
%%     ~~Describe \code{cbs.segs.offspring} here~~
}
  \item{segment.md}{
%%     ~~Describe \code{segment.md} here~~
}
  \item{offspring.segs}{
%%     ~~Describe \code{offspring.segs} here~~
}
  \item{mindistance.threshold}{
%%     ~~Describe \code{mindistance.threshold} here~~
}
  \item{narrowRanges}{
%%     ~~Describe \code{narrowRanges} here~~
}
  \item{calculate.lr}{
%%     ~~Describe \code{calculate.lr} here~~
}
  \item{prOutlier}{
%%     ~~Describe \code{prOutlier} here~~
}
  \item{prMosaic}{
%%     ~~Describe \code{prMosaic} here~~
}
  \item{prob.nonMendelian}{
%%     ~~Describe \code{prob.nonMendelian} here~~
}
  \item{mu.logr}{
%%     ~~Describe \code{mu.logr} here~~
}
  \item{baf.sds}{
%%     ~~Describe \code{baf.sds} here~~
}
  \item{pruneMinimumDistance}{
%%     ~~Describe \code{pruneMinimumDistance} here~~
}
  \item{prune.by.call}{
%%     ~~Describe \code{prune.by.call} here~~
}
  \item{min.coverage}{
%%     ~~Describe \code{min.coverage} here~~
}
  \item{returnEmission}{
%%     ~~Describe \code{returnEmission} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{DNAcopy.verbose}{
%%     ~~Describe \code{DNAcopy.verbose} here~~
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
