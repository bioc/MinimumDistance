\name{callDenovoSegments}
\alias{callDenovoSegments}
\title{
  Posterior calling for segmented data.
}
\description{
  This function provides a convenient wrapper for the segmentation and
  posterior calling steps.
}
\usage{
callDenovoSegments(path = "", filenames, ext = ".txt", pedigreeData, cdfname, chromosome = 1:22, segmentParents, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
    \code{character} string indicating path to BeadStudio files.
}
  \item{filenames}{
    \code{character}: BeadStudio filenames
}
  \item{ext}{
    \code{character}: filename extension
}
  \item{pedigreeData}{
    Object of class \code{Pedigree}.
}
  \item{cdfname}{
    \code{character} indicating the annotation package.
}
  \item{chromosome}{
    \code{integer} indicating which autosomal chromosomes to process
}
\item{segmentParents}{
  \code{Logical}: whether to segment the parental log R ratios using
  circular binary segmentation
}
  \item{verbose}{
    \code{Logical}: whether to display verbose output indicating progress.
}
  \item{\dots}{
    Additional arguments can be specified for the \code{segment}
    function in the package \code{DNAcopy}
}
}

\details{

  A wrapper for the segmentation and posterior calling steps.  Circular
  binary segmentation is performed on the minimum distance and the
  offspring log R ratios.  The minimum distance (the signed minimum of
  the absolute difference of the offspring and parental log2 R ratios)
  should have mean zero in regions of inherited CNV and normal diploid
  genomes.  Hence, the segmentation algorithm should smooth over
  inherited CNV in which the offspring has the same copy number as the
  parent as well as normal diploid regions. Therefore, the posterior
  call for a segment that contains both inherited CNV and normal diploid
  regions will depend on the relative size of these regions.  For
  example, if most of the segment is diploid the posterior call would be
  '333' (see State Symbols below for details regarding the state
  symbols).  Conversely, if most of the region contains a deletion
  transmitted from the mother, the state call will be '322'.  Of primary
  interest are regions for which the offspring copy number differs from
  the parental copy numbers.  Note that for such regions, the CNV in the
  offspring may be Mendelian or non-Mendelian.

  The likelihood for the called state and the likelihood for the normal
  state are named 'lik.state' and 'lik.norm', respectively.  The ratio
  can provide a useful rank whereby large values indicate strong
  evidence for the called state.

  State Symbols:

    The MinimumDistance states can be easily remembered as the latent
  copy number offset by 1 for the father, mother, and offspring,
  respectively. For example, a region of the genome for which all
  individuals in the trio are diploid would have state '333'.  By
  contrast, a region for which the father is diploid and the mother and
  offspring are hemizygous would have the state symbol '322'.

}

\value{ A object of class \code{RangedDataCBS} with the predicted trio
  copy number sates.  Only states for which the offspring copy number
  differs from the parental copy numbers are reliable.  In particular,
  the normal state (state '333') may contain both normal and inherited
  CNV as the minimum distance for both states is near zero.  }


\author{
R. Scharpf
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{

  \code{\link{state}} for accessing the state symbols (posterior calls)
  for the genomic ranges.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
path <- system.file("extdata", package="MinimumDistance")
fnames <- list.files(path, pattern=".txt")
ped <- new("Pedigree",
	   trios=data.frame(F=fnames[1], M=fnames[2], O=fnames[3], stringsAsFactors=FALSE),
	   trioIndex=data.frame(individualId=fnames,
	                        memberId=c("F","M","O"),
                                index.in.pedigree=rep(1,3),
	                        stringsAsFactors=FALSE))
map.segs <- callDenovoSegments(path=path,
			       filenames=fnames,
			       ext="",
			       pedigreeData=ped,
			       cdfname="human610quadv1b",
			       chromosome=1,
			       segmentParents=FALSE)
}
\keyword{htest}
\keyword{models}
\keyword{robust}
\keyword{smooth}
\keyword{classif}

