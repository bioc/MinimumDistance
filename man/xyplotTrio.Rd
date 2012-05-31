\name{xyplotTrio}
\alias{xyplotTrio}

\title{  Function for plotting log R ratios, BAFs, and the minimum distance for
  a trio using lattice.}

\description{
  Function for plotting log R ratios, BAFs, and the minimum distance for
  a trio using lattice.
}
\usage{
xyplotTrio(rd, object, frame = 2e+05, lrr.segments, md.segments, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rd}{ A genomic interval represented as a
    \code{RangedDataCNV}-derived object such as the
    \code{RangedDataHMM} class}

  \item{object}{ An object of class \code{TrioSetList}.}

  \item{frame}{
    The distance (in basepairs) to plot upstream and downstream of the
    genomic interval in \code{rd}.
  }

  \item{lrr.segments}{Optional. A \code{RangedDataCBS} object containing
  the results of a segmentation of the log R ratios (e.g, by circular
  binary segmentation).}

  \item{md.segments}{Optional. A \code{RangedDataCBS} object containing
  the results of a segmentation of the minimum distance (e.g, by circular
  binary segmentation).
}

\item{\dots}{
  Additional arguments are passed to \code{xypanelTrio} and low-level
  lattice functions including \code{panel.xyplot} and \code{lpoints}.}

}


\value{
  A trellis object
}

\author{
R Scharpf
}

\seealso{
  \code{\linkS4class{RangedDataHMM}},
  \code{\linkS4class{RangedDataCBS}}, \code{\link{xyplot}}, \code{\link{xypanelTrio}}
}

\examples{

library(oligoClasses)
data(trioSetListExample)
data(map.segs)
data(lrr.segs)
data(md.segs)
## select a range with a possible de novo copy number alteration in the offspring
rd <- map.segs[which(isDenovo(state(map.segs)))[1], ]
trioSet <- trioSetList[[match(as.character(chromosome(rd)), paste("chr",chromosome(trioSetList),sep=""))]][, match(sampleNames(rd), sampleNames(trioSetList))]
mindist(trioSet) <- calculateMindist(lrr(trioSet))
figs <- xyplotTrio(rd=rd,
		   object=trioSet,
		   frame=200e3,
		   ylab="log R ratio and BAFs",
		   xlab="physical position (Mb)",
		   panel=xypanelTrio,
		   scales=list(x="same",
		   y=list(alternating=1,
		   at=c(-1, 0, log2(3/2), log2(4/2)),
		   labels=expression(-1, 0, log[2](3/2), log[2](4/2)))),
		   lrr.segments=lrr.segs,
		   md.segments=md.segs,
		   layout=c(1, 4),
		   col.hom="grey50",
		   col.het="grey50",
		   col.np="grey20",
		   state.cex=0.8,
		   cex=0.3,
		   ylim=c(-3, 1.5),
		   par.strip.text=list(lines=0.8, cex=0.6),
		   key=list(text=list(c(expression(log[2]("R ratios")), expression("B allele freqencies")),
			    col=c("black", "blue")), columns=2))
print(figs)
}
\keyword{hplot}

