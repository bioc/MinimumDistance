\name{minimumDistancePlot}
\alias{minimumDistancePlot}

\title{
  Creates a list of trellis objects.
}
\description{

  Creates a list of \code{trellis} objects containing information on
  copy number estimates and B allele frequencies.

}
\usage{
minimumDistancePlot(trioSets, ranges, md.segs, cbs.segs, frame = 2e+06, cex = 0.2, scales.cex = 0.5, pch = 21, panelLabelsRight = c("father", "mother", "offspring"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trioSets}{
    A \code{TrioSetList} object.
}
  \item{ranges}{
    A \code{RangedDataCNV} object.
}
  \item{md.segs}{
    A \code{RangedDataCBS} object created from the segmentation of the
    minimum distance.
}
\item{cbs.segs}{

  Optional.  When not missing, \code{cbs.segs} must be a
    \code{RangedDataCBS} object created from the segmentation of the
    copy number estimates for the father, mother, and offspring.

  }

  \item{frame}{Integer.  Indicates how big of a frame around the start
    and end coordinates of a genomic range should be included for
    plotting.  For example, 2 Mb indicates that we will plot copy number
    and B allele frequencies for all markers in the interval [start -
    2Mb, end + 2Mb].}

  \item{cex}{
    Controls size of plotting symbols.
}
  \item{scales.cex}{
    Controls size of axis labels.
}
  \item{pch}{
    Plotting symbol.
}
  \item{panelLabelsRight}{
    Character string for annotating the panels.
}
}
\details{

}
\value{
 A list of trellis objects.
}
\author{
Rob Scharpf
}

\seealso{
\code{\link{gridlayout}}, \code{\link{gridlayout2}}
}
\examples{
	\dontrun{
		data("trioSet-hapmap")
		tlist <- minimumDistancePlot(trioSets=trioSetList,
					     ranges=ranges,
					     md.segs=md.segs,
					     ##cbs.segs=cbs.segs,
					     cex=0.2,
					     scales.cex=0.6)
		## plots of log R ratios
		print(tlist[[1]])
		## plots of B allele frequencies
		print(tlist[[2]])
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{dplot}

