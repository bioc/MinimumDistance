\name{gridlayout}
\alias{gridlayout}
\title{Plot copy number estimates and B allele frequencies}
\description{

  Use lattice graphics to plot the copy number and B allele frequencies
  for a genomic interval.

}
\usage{
gridlayout(filename, lattice.object, rd, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{character string.  Can be missing if nothing is to be saved.}
  \item{lattice.object}{a list of \code{trellis} objects. See details.}
  \item{rd}{A \code{RangedDataCNV} object.  If the object has more than
    1 range, only the first is used.}
  \item{\dots}{Additional arguments to \code{gpar} to control the
    appearance of the vertical lines drawn at the start and stop
    coordinates of the range. See details}
}

\details{

  \code{lattice.object} is a list of trellis objects created for one
  genomic range stored in the \code{rd} object.  The first element in
  the \code{lattice.object} is typically a \code{trellis} object for the
  log R ratios or copy number estimates. The second list element is the
  \code{trellis} object for the B allele frequencies.  Typically the
  \code{lattice.object} includes the data from markers within the
  interval specified by the \code{rd} object, as well as a window that
  frames the region.  See the function \code{\link{minimumDistancePlot}}
  to create the \code{lattice.object}.

  The \code{rd} object is a \code{RangedDataCNV} object for a single
  range.  The code{lattice.object} should have been created from the
  \code{rd} object.

  Vertical lines are drawn at the start and end coordinates in the
  \code{rd} object.  Appearance of the lines can be controlled by
  additional arguments passed to \code{gpar} through \code{\ldots}.

}
\value{
  Nothing is returned.
}
\author{
Rob Scharpf
}

\seealso{

\code{\link{minimumDistancePlot}} for creating the
\code{lattice.object}.

See \code{\link{gpar}} for controlling the appearance of the lines
bordering the range in the \code{rd} object.

}
\examples{
	data(trioSetList)
	data(called.ranges)
	data(md.segs)
	data(cbs.segs)
	cbs.segs$family <- MinimumDistance:::getFamilyName(cbs.segs, trioSetList)
	##trace(minimumDistancePlot, browser)
	##plot the first third range
	trellisList <- minimumDistancePlot(trioSetList=trioSetList,
					   ranges=called.ranges[92,],
					   md.segs=md.segs,
					   cbs.segs=cbs.segs,
					   cex=0.2,
					   scales.cex=0.6)
	logr.trellis.list <- trellisList[[1]]
	baf.trellis.list <- trellisList[[2]]
	gridlayout(lattice.object=list(logr.trellis.list[[1]],
		   baf.trellis.list[[1]]), rd=called.ranges[92,])

	## use a frame of 500kb around the start and end position of the segment
	trellisList <- minimumDistancePlot(trioSetList=trioSetList,
					   ranges=called.ranges[92,],
					   md.segs=md.segs,
					   cbs.segs=cbs.segs,
					   cex=0.2,
					   scales.cex=0.6,
					   frame=500e3)
	logr.trellis.list <- trellisList[[1]]
	baf.trellis.list <- trellisList[[2]]
	gridlayout(lattice.object=list(logr.trellis.list[[1]],
		   baf.trellis.list[[1]]), rd=called.ranges[92,])
}
\keyword{plot}

