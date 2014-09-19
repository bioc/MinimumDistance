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
#' \code{MinDistExperiment} class extends \code{SummarizedExperiment},
#' and so many of the methods defined for \code{SummarizedExperiment}
#' such as \code{findOverlaps} are available through inheritance.
#'
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
