\name{xypanelTrio}
\alias{xypanelTrio}

\title{
  Panel function for a lattice-style plot of log R ratios, BAFs, and the minimum
  distance for a trio.
}

\description{
    Panel function for plotting log R ratios, BAFs, and the minimum
    distance for a trio. This function is not called directly by the
  user, but passes as the argument to panel in the function \code{xyplotTrio}.

}
\usage{
xypanelTrio(x, y, memberId, baf, is.snp, range, lrr.segments,
md.segments, baf.color="blue", col.hom = "grey20", fill.hom = "lightblue", col.het = "grey20", fill.het = "salmon", col.np = "grey20", fill.np = "grey60", show.state = TRUE, cex.state = 1, col.state = "blue", ped, ..., subscripts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Physical position.  Remark: ploted as x/1e6 for Mb units.
}
  \item{y}{
    log R ratios or minimum distance
}
  \item{memberId}{
    father, mother or offspring.
}
  \item{baf}{
    B allele frequencies
}
  \item{is.snp}{
    Logical for whether the marker is at a SNP
  }

  \item{range}{
    A \code{RangedDataCNV}-derived object with physical position for the
    range, a 'chrom', and an 'id' columns.  Should be only one genomic interval.
  }

  \item{lrr.segments}{
    Optional.  A \code{RangedDataCBS} containing the results from a
    (CBS) segmentation of the log R ratios
}
  \item{md.segments}{
    Optional.  A \code{RangedDataCBS} containing the results from a
    (CBS) segmentation of the minimum distance.
  }

  \item{baf.color}{
    Color for plotting BAFs.
  }

  \item{col.hom}{
    Color for homozygous genotypes
}
  \item{fill.hom}{
    Fill color for homozygous genotypes
}
  \item{col.het}{
    Color for heterozygous genotypes
}
  \item{fill.het}{
    Fill color for heterozygous genotypes
}
  \item{col.np}{
    Color for nonpolymorphic markers (genotypes irrelevant)
}
  \item{fill.np}{
    Fill color for nonpolymorphic markers (genotypes irrelevant)
}
  \item{show.state}{
    Whether to display the inferred state from the \code{rd} object.
}
  \item{cex.state}{
    The size of the displayed state. Ignored if \code{show.state} is FALSE.
}
  \item{col.state}{
    Color of the text used for displaying the state.
}

  \item{ped}{
    Object of class \code{Pedigree}
}
  \item{\dots}{
    Additional arguments to \code{lpoints} and \code{panel.xyplot}
}
  \item{subscripts}{
    See \code{xyplot}
}
}

\value{
  This function is used to plot data in the panels of a lattice
    graphic and is not called directly by the user.
}

\author{
R Scharpf
}

\seealso{
  \code{\link{xyplot}}, \code{\link{xyplotTrio}}, \code{\link{lpoints}}, \code{\link{panel.xyplot}}
}
\keyword{dplot}


