\name{read.bsfiles}
\alias{read.bsfiles}

\title{
  Read BeadStudio files
}
\description{
  Reads BeadStudio files  with headers "<ID>.Log.R.Ratio" and "<ID>.B.Allele.Freq"
}
\usage{
read.bsfiles(path = "./", filenames, ext = "", row.names = 1, sep = "\t", as.is = TRUE, header = TRUE, drop = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
    \code{character}: path to BeadStudio files
}
  \item{filenames}{
    \code{character}: filenames
}
  \item{ext}{
    \code{character}: file extension
}
  \item{row.names}{
    \code{integer} indicating which column contains the feature identifiers
}
  \item{sep}{
    passed to \code{read.table}
}
  \item{as.is}{
    passed to \code{read.table}
}
  \item{header}{
    passed to \code{read.table}
}
  \item{drop}{
    \code{logical}:  if TRUE, array is labeled by feature and sample
  identifiers.  If FALSE, \code{dimnames} is NULL.
}
  \item{\dots}{
    Additional arguments passed to \code{read.table}
}
}
%\details{
%
%}
\value{
An array of log R ratios and B allele frequencies
}
%\references{
%%% ~put references to the literature/web site here ~
%}
\author{
R. Scharpf
}
%\note{
%%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{read.table}}
}
\examples{
path <- system.file("extdata", package="MinimumDistance")
fnames <- list.files(path, pattern=".txt")
dat <- read.bsfiles(path=path, fnames, ext="")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}

