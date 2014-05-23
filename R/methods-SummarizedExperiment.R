#' @importMethodsFrom GenomicRanges SummarizedExperiment
computeEmissionProbs <- function(object){
  object <- NA_filter(object)
  emitF <- .emission_one_sample(object[, 1])
  emitM <- .emission_one_sample(object[, 2])
  emitO <- .emission_one_sample(object[, 3])
  cn_states <- paste0("CN", c(0:2, "2-ROH", 3, 4))
  colnames(emitF) <- colnames(emitM) <- colnames(emitO) <- cn_states
  SummarizedExperiment(assays=SimpleList(father=emitF,
                         mother=emitM,
                         offspring=emitO),
                        rowData=rowData(object))
}

.map_snpexp <- function(object,
			 ranges,
                         transition_param,
                         emission_param,
			 mdThr=0.9,...){
  pkgs <- c("GenomicRanges", "VanillaICE", "oligoClasses", "matrixStats", "MinimumDistance")
  build <- genome(object)[1]
  mads <- pmax(ranges$mindist.mad, .1)
  ranges$exceeds.md.thr <- abs(ranges$seg.mean/mads) > mdThr
  ##fit <- hmm2(se) ## A GRangesList
  emitlist <- updateEmission(object)
  browser()
  granges <- sort(ranges)
  ranges <- loglik2(emit=emit,
                    ranges=granges,
                    pr.nonmendelian=pr.nonmendelian,
                    overlapFun=overlapFun)
  chr.arm <- .getArm(chromosome(ranges), start(ranges), build)
  ranges <- combineRangesByFactor(ranges, paste(chr.arm, state(ranges), sep="_"))
  ranges
  results <- unlist(GRangesList(results))
  metadata(results) <- metadata(ranges)
}

setMethod(MAP, c("SnpArrayExperiment", "GRanges"), function(object,
                                                            ranges,
                                                            transition_param=TransitionParam(),
                                                            emission_param=EmissionParam(),
                                                            mdThr=0.9, ...){
  .map_snpexp(object=object,
              ranges=ranges,
              transition_param=transition_param,
              emission_param=emission_param,
              mdThr=mdThr,...)
})
