few_markers <- function(){
  library(data.table)
  library(VanillaICE)
  ##testdir <- "~/Software/bridge/MinimumDistance/inst/unitTests/test_data"

  g <- sort(readRDS(file.path(testdir, "granges.rds")))
  me <- readRDS(file.path(testdir, "me.rds"))
  me <- NA_filter(me)
  g_calls <- readRDS(file.path(testdir, "g_calls.rds"))
  i <- subjectHits(findOverlaps(g_calls[2], me))

  e_param <- EmissionParam(temper=100)
  param <- MinDistParam(thin=20L, emission=e_param)

  b <- baf(me)
  r <- lrr(me)
  assay_list <- list(r[,1], b[,1])
  emissions <- calculateEmission(assay_list, e_param)
  t_param <- TransitionParam(taup=1e10, taumax=1)
  transition_prob <- calculateTransitionProbability(me, t_param)
  hmm_param <- HmmParam(emission=emissions,
                        emission_param=e_param,
                        transition=transition_prob,
                        chrom=chromosome(me),
                        loglik=LogLik(tolerance=0.05))
  while(doUpdate(hmm_param)) hmm_param <- baumWelchUpdate(hmm_param, assay_list)
  i <- subjectHits(findOverlaps(g_calls[5], me))
  le <- log(emission(hmm_param)[i, ])
  colSums(le)
  b[i,]
  r[i,]

  e_param <- EmissionParam(temper=100)
  param <- MinDistParam(thin=20L, emission=e_param)
  load_all("~/Software/bridge/MinimumDistance")
  ##trace(compute_posterior, browser)
  tmp <- compute_posterior(me, g, param, metadata(g)$mad)

  object <- filterExperiment(me, g, param, 0.09)
  object <- NA_filter(object)
  ##computeEmissionProbs(object)

  transition_param <- TransitionParam()
  F <- updateHmmParams(object[, "father"], emission(param),
                       transition_param=transition_param)
  ep <- emissionParam(F)
  e_param <- EmissionParam(cn_means=cn_means(ep), cn_sds=cn_sds(ep),
                           baf_means=baf_means(ep), baf_sds=baf_sds(ep))
  emission(param) <- e_param
  M <- updateHmmParams(object[, "mother"], emission(param), transition_param=transition_param)
  ## Again, update initial values from Mother
  ep <- emissionParam(M)
  e_param <- EmissionParam(cn_means=cn_means(ep), cn_sds=cn_sds(ep),
                           baf_means=baf_means(ep), baf_sds=baf_sds(ep))
  emission(param) <- e_param
  objList <- lapply(offspring(object), function(id, x) x[, id], x=object)
  Olist <- lapply(objList, updateHmmParams, param=emission(param),
                  transition_param=transition_param)

  ##updateHmmParam
  emission_param <- emission(param)
  assay_list <- assayList(object)
  chrom <- as.character(seqnames(object))
  emissions <- calculateEmission(assay_list, emission_param)
  transition_prob <- calculateTransitionProbability(object, transition_param)
  hmm_param <- HmmParam(emission=emissions,
                        emission_param=emission_param,
                        transition=transition_prob,
                        chrom=chrom)
  while(doUpdate(hmm_param)){
    hmm_param <- baumWelchUpdate(hmm_param, assay_list)
  }

  ## baum welch update
  fit <- calculateViterbi(hmm_param)
  LL <- loglik(hmm_param)
  loglik(LL) <- loglik(fit)
  e_param <- updateParam(assay_list, emissionParam(hmm_param), fit)
  emit <- calculateEmission(assay_list, e_param)
  HmmParam(emission=emit,
           emission_param=e_param,
           transition=transition(hmm_param),
           chromosome=hmm_param@chromosome,
           loglik=LL,
           viterbi=fit,
           verbose=verbose(hmm_param))

  ##calculateEmission
  means <- cn_means(param)
  sds <- cn_sds(param)
  limits <- CN_range(param)
  x[[1]] <- threshold(x[[1]], limits)
  emit_cn <- mapply(dnorm, mean=as.list(means), sd=as.list(sds),
                    MoreArgs=list(x=x[[1]]))
  emit_baf <- .calculateBAFEmission(x[[2]], param)
  ## must stay on probability scale for .viterbi2 (do not log)
  ##
  ## Temper the CN emission probabilities because of outliers
  p <- proportionOutlier(param)
  emit_cn <- (1-p)*emit_cn + p*dunif(0, limits[1], limits[2])
  emit_baf <- (1-1e-5)*emit_baf + 1e-5*1
  emit <- emit_cn * emit_baf
  colnames(emit) <- .hmm_states()
}
