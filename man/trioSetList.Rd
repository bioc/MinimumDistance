\name{TrioSetList}
\alias{TrioSetList}

\title{Constructs and object of class \code{TrioSetList}}

\description{Constructs and object of class \code{TrioSetList}}
\usage{
TrioSetList(chromosome=integer(), pedigreeData=Pedigree(), sample.sheet,
row.names=NULL, lrr, baf, featureData, cdfname, ffname="",genome)
}

\arguments{
  \item{chromosome}{Vector of integers indicating which chromosomes are
    stored in the \code{TrioSetList} object.}

  \item{pedigreeData}{A object of class \code{Pedigree}.}

  \item{sample.sheet}{A \code{data.frame} of sample-level covariates.}

  \item{row.names}{The \code{row.names} must match
  the sample identifiers stored in the argument to \code{pedigreeData}.
  If \code{sample.sheet} is missing, \code{row.names} is ignored.}

  \item{lrr}{
    Matrix of log R Ratios
  }

  \item{baf}{ matrix of B allele frequencies}

  \item{featureData}{Object of class \code{AnnotatedDataFrame}
    containing feature annotation. If missing, the argument
    \code{cdfname} must be specified.}

  \item{cdfname}{A character string providing the name of the annotation
    package.}

  \item{ffname}{Prefix for naming ff files. Ignored if \pkg{ff} package
    is not loaded.}

    \item{genome}{character string indicating UCSC genome build. Only
  "hg19" is allowed for annotation packages that support a single
  build. For annotation packages that support multiple builds, valid
  builds are "hg18" and "hg19".}

}

\value{
An object of class \code{TrioSetList}
}

\seealso{
  \code{\linkS4class{TrioSetList}}, \code{\link{phenoData}}
}
\examples{
library(oligoClasses)
TrioSetList()
TrioSetList(chromosome=1:22)

## A more realistic example
## Note that a data.frame containing covariates on the samples can be
## passed through the sample.sheet argument
library(human610quadv1bCrlmm)
path <- system.file("extdata", package="MinimumDistance")
load(file.path(path, "pedigreeInfo.rda"))
load(file.path(path, "sample.sheet.rda"))
load(file.path(path, "logRratio.rda"))
load(file.path(path, "baf.rda"))
nms <- paste("NA",substr(sample.sheet$Sample.Name, 6, 10),sep="")
trioSetList <- TrioSetList(lrr=logRratio, ## must provide row.names
			   baf=baf,
			   pedigree=Pedigree(pedigreeInfo),
			   sample.sheet=sample.sheet,
			   row.names=nms,
			   cdfname="human610quadv1bCrlmm",
			   genome="hg18")
motherPhenoData(trioSetList)
fatherPhenoData(trioSetList)
offspringPhenoData(trioSetList)
## log R ratios for the first trioSetList element
str(lrr(trioSetList)[[1]])
## B allele frequencies for the first trioSetList element
str(baf(trioSetList)[[1]])
}
\keyword{classes}

