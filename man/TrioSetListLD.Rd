\name{TrioSetListLD}
\alias{TrioSetListLD}
\title{
  Constructor for TrioSetList class for large data
}
\description{
  Constructor for TrioSetList class for large data
}
\usage{
TrioSetListLD(path, fnames, ext="", samplesheet, row.names,
pedigreeData, featureData, annotationPkg, outdir=ldPath(), ffprefix="",
genome=c("hg19", "hg18"))
}
\arguments{
  \item{path}{
    Path to plain-text files containing log R ratios and B allele
    frequencies.  Files should contain data for a single sample.
}
  \item{fnames}{
    Character string providing the filenames.
  }

  \item{ext}{
    Character string indicating whether the \code{fnames} has a file
    extension (e.g., ".txt")
  }

  \item{samplesheet}{
    (Optional) \code{data.frame} containing phenotypic / experimental covariates on
    the samples.  Note that if \code{samplesheet} is provided,
    \code{row.names} must be specified.
  }

\item{row.names}{

    Character vector indicating the sample id for each row in
    \code{samplesheet}.  \code{row.names} should be unique and, ideally,
    correspond to \code{fnames},

}
  \item{pedigreeData}{

    An object of class \code{Pedigree}.

  }

  \item{featureData}{ A \code{GenomeAnnotatedDataFrame} }

  \item{annotationPkg}{
    Character string indicating the annotation package used to extract
    information on the features (chromosome, physical position, and
    whether the feature is polymorphic ('isSnp')).
  }

  \item{outdir}{Character string indicating the path for storing \code{ff}
    objects.  Ignored if the \pkg{ff} package is not loaded.}

  \item{ffprefix}{character string indicating the prefix used to name ff
    objects. Ignored if the \pkg{ff} package is not loaded.}

  \item{genome}{character string indicating UCSC genome build. Only
  "hg19" is allowed for annotation packages that support a single
  build. For annotation packages that support multiple builds, valid
  builds are "hg18" and "hg19".}

}
\details{

  If the \pkg{ff} package  is loaded, the \code{assayData} elements will
  be of class \code{ff_array}.  Otherwise, the \code{assayData} elements
  will be ordinary arrays. For large datasets (or for computers with
  limited RAM), loading the \pkg{ff} may be required.
}

\value{
  A \code{TrioSetList} object
}
\author{
R. Scharpf
}

\seealso{
  \code{\linkS4class{TrioSetList}}
}
\examples{
if(require("ff")){
	library(ff)
        library(oligoClasses)
	ldPath(tempdir())
	path <- system.file("extdata", package="MinimumDistance")
	fnames <- list.files(path, pattern="[FMO].txt")
	trioSetListff <- TrioSetListLD(path=path,
				       fnames=fnames,
				       pedigreeData=Pedigree(data.frame(F="F.txt",
				       M="M.txt", O="O.txt")),
				       annotationPkg="human610quadv1bCrlmm",
				       outdir=ldPath(), prefix="test",
				       genome="hg19")
}
}
\keyword{manip}

