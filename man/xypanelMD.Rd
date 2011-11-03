\name{xypanelMD}
\alias{xypanelMD}

\title{Panel function for plotting log R ratios and B allele frequencies
  for a trio}
\description{

  This is a panel function for lattice's xyplot that is useful for
  plotting log R ratios and B allele frequencies for a trio

}
\usage{
xypanelMD(x, y, id, gt, is.snp, range, col.hom = "grey20", fill.hom = "lightblue", col.het = "grey20", fill.het = "salmon", col.np = "grey20", fill.np = "grey60", show.state = TRUE, lrr.segs, md.segs, ..., subscripts)
}

\arguments{
  \item{x}{
    Physical position in MB.
}
  \item{y}{
    Either the log R ratios, minimum distance, or B allele frequency.
}
  \item{id}{
    The id for the trio.
}
  \item{gt}{
    Genotype call (if available)
}
  \item{is.snp}{
    Indicator for whether the marker is polymorphic.
}
  \item{range}{
    A genomic interval of class \code{RangedDataHMM}.
}
  \item{col.hom}{
    Color for homozygous genotypes.
}
  \item{fill.hom}{
    Fill color for homozygous genotypes (ignored for some plotting
  symbols).
}
  \item{col.het}{
    Color for heterozygous genotypes.
}
\item{fill.het}{
    Fill color for heterozygous genotypes (ignored for some plotting
  symbols).

}
  \item{col.np}{
      Color for nonpolymorphic markers
}
\item{fill.np}{

      Fill color for nonpolymorphic markers (ignored for some plotting
  symbols).
}
  \item{show.state}{
    Whether to display the inferred copy number state stored in the
  object passed to the \code{range} argument.
}
  \item{lrr.segs}{
    Object of class \code{RangedDataCBS} containing the results of the
  CBS segmentation of the log R ratios. Optional
}
\item{md.segs}{
      Object of class \code{RangedDataCBS} containing the results of the
  CBS segmentation of the minimum distance. Optional
}
  \item{\dots}{
    Additional arguments to panel.xyplot.
}
  \item{subscripts}{
    See lattice xyplot.
}
}
\details{
  This function is not typically called directly, but is passed to the
  \code{panel} argument in the \code{xyplot} method.
}
\value{
  Nothing is returned.
}

\author{
R. Scharpf
}

\seealso{
  \code{\link{xyplot}}, \code{\link{xypanel}}, \code{\link{panel.xyplot}}
}
\examples{
data(trioSetListExample)
trioSet <- stack(trioSetList)
data(map.segs)
data(lrr.segs)
data(md.segs)
## select a range with a possible de novo copy number alteration in the offspring
rd <- map.segs[which(isDenovo(state(map.segs)))[1], ]

lrrfig <- xyplot(lrr ~ x | id,
		 data=trioSet,
		 range=rd,
		 lrr.segs=lrr.segs,
		 md.segs=md.segs,
		 frame=500e3,
		 panel=xypanelMD,
		 cex=0.2,
		 pch=21,
		 border="orange",
		 scales=list(x="free"),
		 par.strip.text=list(lines=0.8, cex=0.6),
		 col.np="royalblue", fill.np="grey",
		 index.cond=list(4:1),
		 return.data.frame=FALSE)
baffig <- xyplot(baf ~ x | id,
		 data=trioSet,
		 range=rd,
		 lrr.segs=lrr.segs,
		 md.segs=md.segs,
		 frame=500e3,
		 panel=xypanel,
		 cex=0.2,
		 pch=21,
		 border="orange",
		 scales=list(x="free"),
		 par.strip.text=list(lines=0.8, cex=0.6),
		 col.np="royalblue", fill.np="grey",
		 index.cond=list(4:1))
VanillaICE:::arrangeSideBySide(lrrfig, baffig)

}
\keyword{dplot}
\keyword{color}

