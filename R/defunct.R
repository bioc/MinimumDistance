#' Defunct functions/classes/methods in the MinimumDistance package
#'
#' The function, class, or data object you asked is defunct.
#'
#' @aliases coerce,RangedDataCNV,GRanges-class
#' @keywords internal
#' @rdname Defunct
#' @name Defunct
NULL
setAs("RangedDataCNV", "GRanges", function(from, to){
  .Defunct("RangedDataCNV class is defunct")
})





#' Deprecated TrioSetList constructor for large data
#'
#' The TrioSetListLD constructor uses ff objects to handle large
#' datasets.  This function is defunct. Use MinDistExperiment instead.
#'
#' @param path Path to plain-text files containing log R ratios and B
#' allele frequencies.  Files should contain data for a single sample.
#' @param fnames     Character string providing filenames.
#' @param ext Character string indicating whether the \code{fnames}
#' has a file extension (e.g., ".txt")
#' @param samplesheet (Optional) \code{data.frame} containing
#' phenotypic / experimental covariates on the samples.  Note that if
#' \code{samplesheet} is provided, \code{row.names} must be specified.
#' @param row.names Character vector indicating the sample id for each
#' row in \code{samplesheet}.  \code{row.names} should be unique and,
#' ideally, correspond to \code{fnames}
#' @param pedigreeData     An object of class \code{Pedigree}.
#' @param featureData  A \code{GenomeAnnotatedDataFrame}
#' @param annotationPkg Character string indicating the annotation
#' package used to extract information on the features (chromosome,
#' physical position, and whether the feature is polymorphic
#' ('isSnp')).
#' @param outdir Character string indicating the path for storing \code{ff}  objects.  Ignored if the \pkg{ff} package is not loaded.
#' @param ffprefix Character string indicating the prefix used to name
#' ff objects. Ignored if the \pkg{ff} package is not loaded.
#' @param genome character string indicating UCSC genome build. Only
#' "hg19" is allowed for annotation packages that support a single
#' build. Supported builds for most platforms are "hg18" and "hg19".
#' @return A \code{TrioSetList} object
#' @seealso   \code{\linkS4class{TrioSetList}}
#' @export
TrioSetListLD <- function(path, fnames, ext="", samplesheet, row.names,
			  pedigreeData,
			  featureData,
			  annotationPkg, outdir=ldPath(),
			  ffprefix="",
			  genome=c("hg19", "hg18")){
  .Defunct(msg="See MinDistExperiment")
##  if(!is(pedigreeData, "Pedigree")) stop()
##  if(missing(featureData)){
##    fD <- GenomeAnnotatedDataFrameFrom(file.path(path, paste(fnames[1], ext, sep="")), annotationPkg, genome=match.arg(genome))
##    fD <- fD[chromosome(fD) < 23 & !is.na(chromosome(fD)), ]
##  } else {
##    fD <- featureData
##    rm(featureData); gc()
##  }
##  ad <- assayDataListLD(path=path,
##                        pedigree=pedigreeData,
##                        ext=ext,
##                        featureData=fD,
##                        ffprefix=ffprefix)
##  if(!missing(samplesheet)){
##    if(missing(row.names)) stop("if samplesheet is provided, row.names can not be missing.")
##    index <- row.names %in% allNames(pedigreeData)
##    sample.sheet <- samplesheet[index, ]
##    row.names <- row.names[index]
##    offsprPhenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE,
##                                              sample.sheet=sample.sheet,
##                                              which="offspring",
##                                              row.names=row.names)
##    fatherPhenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE,
##                                              sample.sheet=sample.sheet,
##                                              which="father",
##                                              row.names=row.names)
##    motherPhenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE,
##                                              sample.sheet=sample.sheet,
##                                              which="mother",
##                                              row.names=row.names)
##  } else {
##    offsprPhenoData <- annotatedDataFrameFrom(pedigreeData, byrow=FALSE, which="offspring")
##    fatherPhenoData <- annotatedDataFrameFrom(pedigreeData, FALSE, which="father")
##    motherPhenoData <- annotatedDataFrameFrom(pedigreeData, FALSE, which="mother")
##  }
##  uchrom <- unique(chromosome(fD))
##  uchrom <- uchrom[order(uchrom)]
##  featureDataList <- vector("list", length(uchrom))
##  for(i in seq_along(uchrom)) {
##    tmp <- fD[chromosome(fD) == uchrom[i], ]
##    featureDataList[[i]] <- tmp[order(position(tmp)), ]
##  }
##  object <- new("TrioSetList",
##                assayDataList=ad,
##                featureDataList=featureDataList,
##                chromosome=uchrom,
##                pedigree=pedigreeData,
##                fatherPhenoData=fatherPhenoData,
##                motherPhenoData=motherPhenoData,
##                phenoData=offsprPhenoData,
##                genome=match.arg(genome))
##  return(object)
}
