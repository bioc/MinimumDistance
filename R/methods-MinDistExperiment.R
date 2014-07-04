offspringIndex <- function(x) grep("offspring", x)

.calculate_mindist <- function(offspring, father, mother){
  d1 <- offspring - father
  d2 <- offspring - mother ##offspring - mother
  I <- as.numeric(abs(d1) <= abs(d2))
  I*d1 + (1-I)*d2
}

.mindist <- function(object){
  ## assumes FMO ordering
  cn <- object$data[["cn"]]
  F <- cn[, "father"]
  M <- cn[, "mother"]
  O <- cn[, offspringIndex(colnames(cn)), drop=FALSE]
  apply(O, 2, .calculate_mindist, father=F, mother=M)
}

.setColnames <- function(object=nm, nm){
  colnames(object) <- nm
  object
}

.mindistnames <- function(x) paste0("mindist_", x[offspringIndex(x)])

#' @importMethodsFrom BiocGenerics colnames
setMethod("colnames", "ShallowSimpleListAssays",
          function (x, do.NULL = TRUE, prefix = "col") {
            colnames(x$data[["cn"]])
          })

setGeneric("MinDistExperiment", function(object, rowData, cn, baf, colData, ...)
           standardGeneric("MinDistExperiment"))

## do nothing
setMethod(SnpGRanges, "SnpGRanges", function(object, isSnp) return(object))

.constructMDE <- function(assays, rowData, colData){
   md <- .setColnames(.mindist(assays), .mindistnames(colnames(assays)))
   se <- new("MinDistExperiment", assays=assays,
             rowData=SnpGRanges(rowData), colData=colData,
             mindist=md)
}

setMethod("MinDistExperiment", c("missing", "GRanges", "matrix", "matrix"),
          function(object, rowData, cn, baf, colData){
            filenames <- colData$filename
            if(!all(colnames(cn) %in% filenames)) stop("colnames of cn matrix are not in the pedigree vector.")
            cn <- cn[, filenames]
            baf <- baf[, filenames]
            colnames(cn) <- colnames(baf) <- rownames(colData)
            assays <- snpArrayAssays(cn=cn, baf=baf)
            .constructMDE(assays, rowData, colData)
          })

#' @export
setMethod("MinDistExperiment", c("ShallowSimpleListAssays", "GRanges"),
          function(object, rowData, cn, baf, colData){
            .constructMDE(object, rowData, colData)
          })

setMethod("MinDistExperiment", c("FileViews", "missing"),
          function(object, rowData, cn, baf, colData){
            al <- assays(object)
            extdata <- system.file("extdata", package=object@annot_pkg)
            load(file.path(extdata, paste0("cnProbes_", build, ".rda")))
            load(file.path(extdata, paste0("snpProbes_", build, ".rda")))
            annot <- rbind(snpProbes, cnProbes)
            fid <- rownames(al$data[["cn"]])
            annot <- annot[fid, ]
            is_snp <- fid %in% rownames(snpProbes)
            rowdata <- GRanges(paste0("chr", annot[, "chr"]), IRanges(annot[, "position"], width=1))
            coldata <- setNames(DataFrame(filenames(object)), "filename")
            .constructMDE(al, rowData=SnpGRanges(rowdata, is_snp), coldata)
          })


setMethod("show", "MinDistExperiment", function(object){
  callNextMethod(object)
  ##cat("MAD(minimum distance): ", round(mad(mindist(object),na.rm=TRUE),2),  "\n")
})

##setMethod("pedigree", "MinDistExperiment", function(object) object$pedigree)

setMethod("mindist", "MinDistExperiment", function(object) object@mindist)

setValidity("MinDistExperiment", function(object){
  msg <- TRUE
  if(!"filename" %in% names(colData(object))){
    msg <- "filename must be in colData"
    return(msg)
  }
  msg
})

setMethod("[", "MinDistExperiment", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    if(is(i, "Rle")) i <- as.logical(i)
    if(is(i, "character")) i <- match(i, rownames(x))
    x@mindist <- x@mindist[i, , drop=FALSE]
  }
  callNextMethod(x, i, j, ..., drop=drop)
})

setMethod("offspring", "MinDistExperiment", function(object) colnames(object)[offspringIndex(colnames(object))])
##setMethod("father", "MinDistExperiment", function(object) colnames(object)["father"])
##setMethod("mother", "MinDistExperiment", function(object) pedigree(object)["mother"])


setMethod("subsetAndSort", "MinDistExperiment",
          function(object, autosomes=seqlevels(object)[1:22]){
            object <- object[chromosome(object) %in% autosomes, ]
            seqlevels(rowData(object), force=TRUE) <- autosomes
            object <- sort(object)
            object <- removeDuplicateMapLoc(object)
            object
          })

removeDuplicateMapLoc <- function(object){
  chr_pos <- paste(chromosome(object), start(object), sep="_")
  is_dup <- duplicated(chr_pos)
  if(any(is_dup)) object <- object[!is_dup, ]
  object
}

##.emission_one_sample <- function(object, param=MinDistParam()){
##  hmm_param <- updateHmmParams(object, emisson_param=emission(param),
##                               transition_param=TransitionParam())
##}
  ##obj <- NA_filter(object)
  ##t_param <- TransitionParam()
  ##transition_probs <- calculateTransitionProbability(obj, t_param)
##  hmm_param <- updateHmmParams(object, emisson_param=emission(param),
##                               transition_prob)
##  if(ncol(object) > 1) stop()
##  r <- drop(lrr(obj))
##  b <- drop(baf(obj))
##  names(r) <- names(b) <- NULL
##  rb <- list(r, b)
##  ##
##  ## should we expose these parameters?
##  ##t_param <- TransitionParam(taup=1e10, taumax=1)
##
##  e_param <- emission(param)
##  ##
##  ##
##  emissions <- calculateEmission(rb, e_param)
##
##  hmm_param <- HmmParam(emission=emissions,
##                        emission_param=e_param,
##                        transition=transition_prob,
##                        chrom=chromosome(object))
##  while(doUpdate(hmm_param)){
##    hmm_param <- baumWelchUpdate(hmm_param, rb)
##  }
##  hmm_param
##}


# #' @importMethodsFrom GenomicRanges SummarizedExperiment
computeEmissionProbs <- function(object, param=MinDistParam()){
  object <- NA_filter(object)
  transition_param <- TransitionParam()
  F <- updateHmmParams(object[, "father"], emission(param), transition_param=transition_param)
  M <- updateHmmParams(object[, "mother"], emission(param), transition_param=transition_param)
  objList <- lapply(offspring(object), function(id, x) x[, id], x=object)
  Olist <- lapply(objList, updateHmmParams, param=emission(param), transition_param=transition_param)
  emitO <- lapply(Olist, emission)
  tmp <- SimpleList(father=emission(F), mother=emission(M))
  tmp2 <- SimpleList(emitO)
  tmp2@listData <- setNames(tmp2@listData, offspring(object))
  SummarizedExperiment(assays=c(tmp, tmp2), rowData=rowData(object))
}



setMethod(MAP2, c("MinDistExperiment", "MinDistGRanges"), function(object, mdgr, param, ...){
  emissions_SE <- computeEmissionProbs(object)
  grl <- mindist(mdgr)  ## will be a GRangesList with length equal to the number of offspring
  mads <- mad(mdgr)[names(grl)]
  grl <- foreach(g=grl, md.mad=mads) %do% {
    compute_loglik(emissions_SE, md_gr=g, param=param, md.mad=md.mad)
  }
  GRangesList(grl)
})

setMethod(MAP2, c("MinDistExperiment", "GRanges"), function(object, mdgr, param, ...){
  object2 <- computeEmissionProbs(object)
  md.mad <- mad(as.numeric(mindist(object)), na.rm=TRUE)
  g <- compute_loglik(object2, md_gr=mdgr, param=param, md.mad=md.mad)
})

segMeanAboveThr <- function(mean, mad, nmad) abs(mean)/max(mad, 0.1) > nmad

posteriorSummaries <- function(log_prior.lik){
  reference <- log_prior.lik[["222"]]
  probs <- exp(log_prior.lik)
  probs <- probs[!is.na(probs)]
  loglik <- log_prior.lik[which.max(log_prior.lik)]
  posterior_log_odds <- posteriorLogOdds(probs)
  ##results$loglik[i] <- round(loglik, 2)
  x <- list(call=names(loglik),
            posterior_log_RR=round(loglik - reference, 2),
            posterior_log_odds=round(posterior_log_odds, 2),
            log_posterior_MAP=round(loglik, 2),
            log_posterior_222=round(log(exp(reference)/(sum(probs))), 2))
  return(x)
}



compute_loglik <- function(object, md_gr, param, md.mad){
  pennparam <- penncnv(param)
  hits <- findOverlaps(md_gr, rowData(object))
  feature_index <- subjectHits(hits)
  results <- .data_frame_posteriorSummaries(md_gr)
  above_thr <- segMeanAboveThr(mean=md_gr$seg.mean, mad=md.mad, nmad=nMAD(param))
  log_emit <- emissionArray(object, epsilon=minimum_emission(pennparam))
  which.offspring <- which(names(assays(object))==gsub("mindist_", "", md_gr$sample[1]))
  log_emit <- log_emit[, c(1, 2, which.offspring), ]
  trio_states <- state(pennparam)
  state.prev <- NULL
  states <- stateNames(pennparam)
  qhits <- queryHits(hits)
  nmarkers <- table(qhits)
  Index <- intersect(as.integer(names(nmarkers)), which(above_thr))
  for(k in seq_along(Index)){
    i <- Index[k]
    index <- feature_index[qhits == i]
    ##  likelihood of the observed data given the model
    LLT <- cumulativeLogLik(log_emit[index, , , drop=FALSE])
    ##  log(likelihood * prior models):
    log_prior.lik <- setNames(sapply(states, posterior,
                                     param=pennparam,
                                     state.prev=state.prev,
                                     log.lik=LLT), states)
    post <- posteriorSummaries(log_prior.lik)
    results$call[i] <- post[["call"]]
    results$log_posterior_MAP[i] <- post[["log_posterior_MAP"]]
    results$log_posterior_222[i] <- post[["log_posterior_222"]]
    results$posterior_log_RR[i] <- post[["posterior_log_RR"]]
    results$posterior_log_odds[i] <- post[["posterior_log_odds"]]
    state.prev <- post[["call"]]
  }
  mcols(md_gr) <- cbind(mcols(md_gr), results)
  md_gr[!is.na(md_gr$call)]
}

posteriorLogOdds <- function(x){
  ## posterior odds that state is denovo
  p <- sum(x[c("220", "221")])/sum(x)
  log(p) - log(1-p)
}



## how to you export a coercion from TrioSet to SnpArrayExperiment?
##  export
##setAs("TrioSet", "SnpArrayExperiment", function(from, to){
##  ped <- pedigree(from)
##  cn <- lrr(from)[, 1, ]/100
##  b <- baf(from)[, 1, ]/1000
##  colnames(b) <- colnames(cn) <- trios(ped)[1, ]
##  gd <- GRanges(paste0("chr", chromosome(from)), IRanges(position(from),
##                                                         width=1),
##                isSnp=isSnp(from))
##  rowdata <- SnpGRanges(gd)
##  se <- SnpArrayExperiment(cn=cn, baf=b, rowData=rowdata)
##  se
##})


setAs("TrioSetList", "MinDistExperiment", function(from, to){
  trioSet <- stack(from)
  as(trioSet, "MinDistExperiment")
})

setAs("TrioSet", "MinDistExperiment", function(from, to){
  if(ncol(from) > 1) warning("only coercing first trio in TrioSet to MinDistExperiment")
  ##trioSet <- stack(trioSetList)[, 1]
  from <- from[, 1]
  ped <- trios(pedigree(from))
  trios <- setNames(as.character(ped), c("father", "mother", "offspring"))
  gd <- GRanges(paste0("chr", chromosome(from)),
                IRanges(position(from),
                        width=1),
                isSnp=isSnp(from))
  r <- .setColnames(lrr(from)[, 1, ], trios)/100
  b <- .setColnames(baf(from)[, 1, ], trios)/1000
  me <- MinDistExperiment(cn=r,
                          baf=b,
                          rowData=gd,
                          colData=setNames(DataFrame(trios), "filename"))
  me
})
