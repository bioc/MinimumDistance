#' De novo copy number alterations in parent-offspring trios
#'
#' @docType package
#' @name MinimumDistance
#' @import ff
#' @import Biobase
#' @importFrom lattice xyplot
#' @importFrom foreach foreach %do% %dopar%
#' @importMethodsFrom oligoClasses genomeBuild
#' @importMethodsFrom oligoClasses state lrr "lrr<-" baf "baf<-" position checkOrder isSnp chromosome GenomeAnnotatedDataFrameFrom
#' @importMethodsFrom IRanges  nrow ncol unlist split dim width length order elementLengths start end
#' @importMethodsFrom IRanges  findOverlaps subjectHits queryHits disjoin countOverlaps
#' @importClassesFrom IRanges DataFrame SimpleList List
#' @importClassesFrom S4Vectors DataTable DataTableORNULL Vector Annotated
#' @importMethodsFrom S4Vectors mcols "mcols<-" values metadata "metadata<-"
#' @importClassesFrom GenomicRanges SummarizedExperiment
#' @importMethodsFrom GenomicRanges rowData colData assays "rowData<-"
#' @import VanillaICE
#' @importFrom GenomicRanges GRangesList GRanges
NULL

# @importClassesFrom Biobase eSet Versioned VersionedBiobase
# @impotMethodsFrom Biobase featureNames
#  @import Biobase
#  @import VanillaICE
#  @importClassesFrom VanillaICE SnpArrayExperiment SnpGRanges
#  @importMethodsFrom VanillaICE SnpGRanges SnpArrayExperiment NA_filter NA_index calculateEmission calculateTransitionProbability
