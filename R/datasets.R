#' An example \code{TrioSetList} object
#'
#' A dataset containing log R ratios and B allele frequencies for two
#' chromosomes, organized as a \code{TrioSetList}.  Each element in
#' the list class is a \code{TrioSet} object.  Both \code{TrioSetList}
#' and \code{TrioSet} classes are deprecated; the example data will be
#' removed in a future release.
#'
#' @format a \code{TrioSetList}
#' @aliases trioSetList
#' @name exampleTrioSetList
NULL

#' An example \code{MinDistExperiment}
#'
#' This dataset contains log R ratios and B allele frequencies from a
#' parent-offspring trio (three individuals). Only markers from
#' chromosomes 7 and 22 are included in this object. The
#' \code{MinDistExperiment} class extends \code{RangedSummarizedExperiment},
#' and so many of the methods defined for \code{RangedSummarizedExperiment}
#' such as \code{findOverlaps} are available through inheritance.
#' @examples
#' \dontrun{
#'     library(oligoClasses)
#'     library(VanillaICE)
#'     library(data.table)
#'     library(BSgenome.Hsapiens.UCSC.hg18)
#'     extdir <- system.file("extdata", package="VanillaICE")
#'     features <- suppressWarnings(fread(file.path(extdir, "SNP_info.csv")))
#'     fgr <- GRanges(paste0("chr", features$Chr), IRanges(features$Position, width=1),
#'                    isSnp=features[["Intensity Only"]]==0)
#'     fgr <- SnpGRanges(fgr)
#'     names(fgr) <- features[["Name"]]
#'     sl <- seqlevels(BSgenome.Hsapiens.UCSC.hg18)
#'     seqlevels(fgr) <- sl[sl %in% seqlevels(fgr)]
#'     seqinfo(fgr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg18)[seqlevels(fgr),]
#'     fgr <- sort(fgr)
#'     files <- list.files(extdir, full.names=TRUE, recursive=TRUE, pattern="FinalReport")
#'     ## parse files
#'     parsedDir <- "ParsedFiles"
#'     if(!file.exists(parsedDir)) dir.create(parsedDir)
#'     views <- ArrayViews(rowRanges=fgr, sourcePaths=files, parsedPath=parsedDir)
#'     dat <- fread(files[1])
#'     select_columns <- match(c("SNP Name", "Allele1 - AB", "Allele2 - AB",
#'                               "Log R Ratio", "B Allele Freq"), names(dat))
#'     index_genome <- match(names(fgr), dat[["SNP Name"]])
#'     scan_params <- CopyNumScanParams(index_genome=index_genome, select=select_columns,
#'                                      cnvar="Log R Ratio",
#'                                      bafvar="B Allele Freq",
#'                                      gtvar=c("Allele1 - AB", "Allele2 - AB"))
#'     invisible(sapply(views, parseSourceFile, param=scan_params))
#'     ped_hapmap <- ParentOffspring(id = "hapmap", father="12287_03",
#'                                   mother="12287_02",
#'                                   offspring="12287_01",
#'                                   parsedPath=parsedPath(views))
#'     ped_list <- ParentOffspringList(pedigrees=list(
#'                                       ParentOffspring(id = "hapmap", father="12287_03",
#'                                                       mother="12287_02",
#'                                                       offspring="12287_01",
#'                                                       parsedPath=parsedPath(views)),
#'                                       ParentOffspring(id = "cleft",
#'                                                       father="22169_03",
#'                                                       mother="22169_02",
#'                                                       offspring="22169_01",
#'                                                       parsedPath=parsedPath(views))))
#'     sample_info <- read.csv(file.path(extdir, "sample_data.csv"), stringsAsFactors=FALSE)
#'     ind_id <- setNames(gsub(" ", "", sample_info$IndividualID), sample_info$File)
#'     colnames(views) <- ind_id[gsub(".csv", "", colnames(views))]
#'     md_exp <- MinDistExperiment(views, pedigree=ped_list[[2]])
#'     seqlevels(md_exp, force=TRUE) <- "chr22"
#'     params <- MinDistParam()
#'     md_gr <- segment2(md_exp, params)
#'     save(md_exp, file="~/Software/bridge/MinimumDistance/data/md_exp.rda")
#'     save(md_gr, file="~/Software/bridge/MinimumDistance/data/md_gr.rda")
#' }
#' @format a \code{MinDistExperiment}
#' @name md_exp
NULL

#' An example \code{MinDistGRanges} object
#'
#' Prior to inferring de novo trio copy number states, the log R
#' ratios are segmented independently for each individual in a
#' \code{ParentOffsping} class. The segmentation results are recorded
#' in separate \code{GRanges} objects for the parents. For
#' segmentation of the offspring log R ratios and the minimum
#' distance, the segments are stored in separate \code{GRangesList}
#' objects.  For convenience, these \code{GRanges},
#' \code{GRangesList}, and pedigree information are bound in a single
#' container referred to as a \code{MinDistGRanges} object.  The
#' example \code{MinDistGRanges} object provided in this package was
#' obtained from the segmentation of the data stored in the example
#' \code{MinDistExperiment} object.
#'
#' @format a \code{MinDistGRanges} object
#' @name md_gr
NULL
