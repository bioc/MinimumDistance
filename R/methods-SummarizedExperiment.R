.emission_one_sample <- function(object){
  if(ncol(object) > 1) stop()
  obj <- NA_filter(object)
  r <- drop(lrr(obj))
  b <- drop(baf(obj))
  e_param <- EmissionParam()
  rb <- list(r, b)
  emissions <- calculateEmission(rb, e_param)
  t_param <- TransitionParam(taup=1e10, taumax=1)
  transition_prob <- calculateTransitionProbability(obj, t_param)
  hmm_param<- HmmParam(emission=emissions, transition=transition_prob)
  LL <- rep(NA, 10)
  delta <- 1
  i <- 1
  while(delta > 1){
    fit <- viterbi(hmm_param)
    LL[i] <- loglik(fit)
    e_param <- updateParam(rb, e_param, fit)
    emission(hmm_param) <- calculateEmission(rb, e_param)
    if(i > 1) delta <- LL[i]-LL[i-1]
    if(i == 10) {
      warning("Log lik ratio never less than 0.05")
      break()
    }
    i <- i+1
  }
  emission(hmm_param)
}


# #' @importMethodsFrom GenomicRanges SummarizedExperiment
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
