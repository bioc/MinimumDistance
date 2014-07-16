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
  F <- cn[, 1]
  M <- cn[, 2]
  ##O <- cn[, offspringIndex(colnames(cn)), drop=FALSE]
  O <- cn[, -(1:2), drop=FALSE]
  apply(O, 2, .calculate_mindist, father=F, mother=M)
}

.setColnames <- function(object=nm, nm){
  colnames(object) <- nm
  object
}

##.mindistnames <- function(x) paste0("mindist_", x[offspringIndex(x)])


setMethod("colnames", "ShallowSimpleListAssays",
          function (x, do.NULL = TRUE, prefix = "col") {
            colnames(x$data[["cn"]])
          })

## do nothing
setMethod(SnpGRanges, "SnpGRanges", function(object, isSnp) return(object))

.set_md_names <- function(id) paste0("md_", id)
.get_md_names <- function(object) .set_md_names(offspring(object))

.constructMDE <- function(assays, rowData, colData, pedigree){
  md <- .setColnames(.mindist(assays), .set_md_names(offspring(pedigree)))
  new("MinDistExperiment",
      assays=assays,
      rowData=SnpGRanges(rowData),
      colData=colData,
      mindist=md,
      pedigree=pedigree)
}

##setMethod("MinDistExperiment", c("missing", "GRanges", "matrix", "matrix"),
##          function(object, rowData, cn, baf, colData){
##            filenames <- colData$filename
##            if(!all(colnames(cn) %in% filenames)) stop("colnames of cn matrix are not in the pedigree vector.")
##            cn <- cn[, filenames]
##            baf <- baf[, filenames]
##            colnames(cn) <- colnames(baf) <- rownames(colData)
##            assays <- snpArrayAssays(cn=cn, baf=baf)
##            .constructMDE(assays, rowData, colData)
##          })
##
##setMethod("MinDistExperiment", c("missing", "missing", "missing", "missing"),
##          function(object, rowData, cn, baf, colData){
##            new("MinDistExperiment")
##          })
##
###' @export
##setMethod("MinDistExperiment", c("ShallowSimpleListAssays", "GRanges"),
##          function(object, rowData, cn, baf, colData){
##            .constructMDE(object, rowData, colData)
##          })
##
##setMethod("MinDistExperiment", c("FileViews", "missing"),
##          function(object, rowData, cn, baf, colData){
##            al <- assays(object)
##            extdata <- system.file("extdata", package=object@annot_pkg)
##            load(file.path(extdata, paste0("cnProbes_", build, ".rda")))
##            load(file.path(extdata, paste0("snpProbes_", build, ".rda")))
##            annot <- rbind(snpProbes, cnProbes)
##            fid <- rownames(al$data[["cn"]])
##            annot <- annot[fid, ]
##            is_snp <- fid %in% rownames(snpProbes)
##            rowdata <- GRanges(paste0("chr", annot[, "chr"]), IRanges(annot[, "position"], width=1))
##            coldata <- setNames(DataFrame(filenames(object)), "filename")
##            .constructMDE(al, rowData=SnpGRanges(rowdata, is_snp), coldata)
##          })


setMethod("MinDistExperiment", c("ArrayViews", "ParentOffspring"),
          function(object=ArrayViews(),
                   pedigree=ParentOffspring(), ...){
            object <- dropSexChrom(object)
            object <- sort(object)
            object <- dropDuplicatedMapLocs(object)
            al <- assays(object)
            .constructMDE(al, rowData=SnpGRanges(rowData(object)), colData=colData(object),
                          pedigree=pedigree)
          })

setMethod("assays", "ArrayViews", function(x, ...){
  ##if(nrow(x) > 1) stop("object contains several trios. Use [i, ] to subset the ith trio prior to calling assays")
  r <- lrr(x)
  b <- baf(x)
  colnames(r) <- colnames(b) <- colnames(x)
  snpArrayAssays(cn=r, baf=b)
})

setMethod("show", "MinDistExperiment", function(object){
  callNextMethod(object)
  ##cat("MAD(minimum distance): ", round(mad(mindist(object),na.rm=TRUE),2),  "\n")
})

setMethod("pedigree", "MinDistExperiment", function(object) object@pedigree)
setReplaceMethod("pedigree", "MinDistExperiment", function(object,value) {
  object@pedigree <- value
  object
})

setMethod("mindist", "MinDistExperiment", function(object) object@mindist)

setReplaceMethod("mindist", "MinDistExperiment", function(object, value) {
  object@mindist <- value
  object
})

##setValidity("MinDistExperiment", function(object){
##  msg <- TRUE
##  if(!"filename" %in% names(colData(object))){
##    msg <- "filename must be in colData"
##    return(msg)
##  }
##  msg
##})

setMethod("[", "MinDistExperiment", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    if(is(i, "Rle")) i <- as.logical(i)
    if(is(i, "character")) i <- match(i, rownames(x))
    x@mindist <- x@mindist[i, , drop=FALSE]
  }
  if(!missing(j)){
    ## special operations when selecting offspring
    ids <- colnames(x)[j]
    if(any(ids %in% offspring(x))){
      offspr <- offspring(x)[offspring(x) %in% ids]
      ped <- pedigree(x)
      ped@offspring <- offspr
      pedigree(x) <- ped
      mindist(x) <- mindist(x)[, .set_md_names(offspr), drop=FALSE]
    }
  }
  callNextMethod(x, i, j, ..., drop=drop)
})

setMethod("offspring", "MinDistExperiment", function(object) offspring(pedigree(object)))


setMethod("father", "MinDistExperiment", function(object) father(pedigree(object)))
setMethod("mother", "MinDistExperiment", function(object) mother(pedigree(object)))


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


computeEmissionProbs <- function(object, param=MinDistParam()){
  ##object <- NA_filter(object)
  transition_param <- TransitionParam()
  F <- updateHmmParams(object[, father(object)], emission(param), transition_param=transition_param)
  ## use emission parameters (updated by Baum Welch) as initial values for Mother
  emission(param) <- emissionParam(F)
  M <- updateHmmParams(object[, mother(object)], emission(param), transition_param=transition_param)
  ## Again, update initial values from Mother
  emission(param) <- emissionParam(M)
  Olist <- list()
  offspr <- offspring(object)
  for(j in seq_along(offspr)){
    id <- offspr[j]
    Olist[[j]] <- updateHmmParams(object[, id], emission(param), transition_param=transition_param)
  }
  emitO <- lapply(Olist, emission)
  tmp <- SimpleList(father=emission(F), mother=emission(M))
  tmp2 <- SimpleList(emitO)
  tmp2@listData <- setNames(tmp2@listData, offspring(object))
  tmp@listData <- setNames(tmp@listData, c(father(object), mother(object)))
  SummarizedExperiment(assays=c(tmp, tmp2), rowData=rowData(object))
}

setGeneric("colMads", function(x, centers=colMedians(x, ...), constant=1.4826, ...)
           standardGeneric("colMads"))


setMethod("colMads", signature(x="MinDistExperiment"),
          function(x, centers=colMedians(x, ...), constant=1.4826, ...){
            colMads(mindist(x), na.rm=TRUE)
          })

segMeanAboveThr <- function(mean, mad, nmad) abs(mean)/max(mad, 0.1) > nmad

posteriorSummaries <- function(log_prior.lik){
  ## avoid numerical issues
  loglik <- log_prior.lik-max(log_prior.lik, na.rm=TRUE)
  map <- loglik[which.max(loglik)]## or which one is zero
  ##range(log_prior.lik)
  reference <- loglik[["222"]]
  probs <- exp(loglik)/sum(exp(loglik), na.rm=TRUE)
  probs <- probs[!is.na(probs)]
  pr222 <- probs["222"]
  pr221 <- probs["221"]
  prMAP <- probs[names(map)]
  p <- sum(probs[c("220", "221")])/sum(probs)
  posterior_log_odds <- log(p) - log(1-p)
  x <- list(call=names(map),
            posteriors=c(round(log(prMAP/pr222), 3),
              round(posterior_log_odds, 3),
              round(prMAP, 3),
              round(pr222, 3),
              round(pr221, 3)))
  return(x)
}

setMethod(MAP2, c("MinDistExperiment", "MinDistGRanges"), function(object, mdgr, param, ...){
  obj <- computePosterior(object, granges=mindist(mdgr), param=param)
  GRangesList(obj)
})

setMethod(MAP2, c("MinDistExperiment", "GRangesList"), function(object, mdgr, param, ...){
  obj <- computePosterior(object, granges=mdgr, param=param)
  GRangesList(obj)
})

setMethod(MAP2, c("MinDistExperiment", "GRanges"), function(object, mdgr, param, ...){
  obj <- computePosterior(object, granges=mdgr, param=param)
  GRangesList(obj)
})


filterIndexForGRanges <- function(object, granges, param){
  mads <- setNames(colMads(object), .get_md_names(object))
  m <- mads[granges$sample[1]]
  above_thr <- segMeanAboveThr(mean=granges$seg.mean, mad=m, nmad=nMAD(param))
  which(above_thr)
}

.compute_trio_posterior <- function(object, granges, param, lemit){
  states <- stateNames(param)
  Index <- filterIndexForGRanges(object, granges, param)
  hits <- findOverlaps(granges, object)
  qhits <- queryHits(hits)
  feature_index <- subjectHits(hits)
  chrom <- chromosome(granges)
  POST <- as.matrix(.mcols_posteriors(granges[Index]))
  calls <- rep(NA, length(Index))
  state.prev <- NULL
  ## This is a nested for loop.  Port to C.
  for(k in seq_along(Index)){
    i <- Index[k]
    index <- feature_index[qhits == i]
    ##  likelihood of the observed data given the model
    LLT <- cumulativeLogLik(lemit[index, , , drop=FALSE])
    ##  log(likelihood * prior models):
    log_prior.lik <- setNames(sapply(states, posterior,
                                     param=penncnv(param),
                                     state.prev=state.prev,
                                     log.lik=LLT), states)
    tmp <- posteriorSummaries(log_prior.lik)
    POST[k, ] <- tmp$posteriors
    calls[k] <- tmp$call
    state.prev <- previousState(Index, k, i, tmp$call, chrom)
  }
  post <- DataFrame(POST)
  granges <- granges[Index]
  mcols(granges) <- c(mcols(granges), post)
  granges$calls <- calls
  mg <- as(granges, "MDRanges")
  versions <- c(packageVersion("VanillaICE"),
                packageVersion("MinimumDistance"))
  metadata(mg) <- list(versions=setNames(versions, c("VanillaICE", "MinimumDistance")))
  mg
}

##setMethod("computePosterior", c("MinDistExperiment", "GRanges"), function(object, granges, param){
##  if(nrow(object) == 0) return(MDRanges())
##  emissions_object <- computeEmissionProbs(object, param)
##  log_emit <- logEmissionArray(emissions_object)
##
##}

setMethod("computePosterior", c("MinDistExperiment", "GRanges"), function(object, granges, param){
  if(nrow(object) == 0) return(MDRanges())
  ns <- unique(granges$sample)
  if(length(ns) > 1) stop("granges$sample must be unique.")
  if(ncol(object) > 3){
    j <- match(gsub("md_", "", granges$sample[1]), colnames(object))
    object <- object[, c(1,2,j)]
  }
  emissions_object <- computeEmissionProbs(object, param)
  log_emit <- logEmissionArray(emissions_object)
  .compute_trio_posterior(object=object, granges=granges,
                          param=param, lemit=log_emit)
})

setMethod("computePosterior", c("MinDistExperiment", "GRangesList"), function(object, granges, param){
  if(nrow(object) == 0) return(MDRanges())
  ## compute emission probabilities
  emissions_object <- computeEmissionProbs(object, param)
  log_emit <- logEmissionArray(emissions_object)
  offspr <- offspring(object)
  J <- foreach(id=offspr) %do% c(1:2, match(id, colnames(object)))
  md_rangesList <- foreach(j = J, g=granges) %do% {
    .compute_trio_posterior(object=object[, j],
                            granges=g,
                            param=param,
                            lemit=log_emit[, j, ])
  }
  md_rangesList
 })




 previousState <- function(Index, k, i, call, chrom){
   if(i == max(Index)) return(call)
   state.prev <- call
   i.next <- Index[k+1]
   if(i.next - i > 1 || chrom[i] != chrom[i.next]) state.prev <- NULL
   state.prev
 }



 ##setMethod("posteriorLogOdds", "numeric", function(object){
 ##  ## posterior odds that state is denovo
 ##  p <- sum(object[c("220", "221")])/sum(object)
 ##  log(p) - log(1-p)
 ##})



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



setMethod("filterExperiment", c("MinDistExperiment", "GRanges"),
          function(object, granges, param){
            .filter_mdexperiment(object, granges, param)
          })

setMethod("filterExperiment", c("MinDistExperiment", "GRangesList"),
          function(object, granges, param){
            g <- unlist(granges)
            .filter_mdexperiment(object, g, param)
          })

setMethod("filterExperiment", c("MinDistExperiment", "MinDistGRanges"),
          function(object, granges, param){
            g <- unlist(mindist(granges))
            .filter_mdexperiment(object, g, param)
          })


.filter_mdexperiment <- function(object, granges, param){
  object <- NA_filter(object)
  mads <- setNames(colMads(object), .get_md_names(object))
  mads <- mads[granges$sample]
  ##mads <- mad(mdgr)[names(mindist(mdgr))]
  above_thr <- segMeanAboveThr(mean=granges$seg.mean, mad=mads, nmad=nMAD(param))
  ##g <- mindist(mdgr)[[1]]
  index <- unique(subjectHits(findOverlaps(granges[above_thr], rowData(object), maxgap=50e3)))
  if(length(index) > 0){
    index2 <- seq_along(rowData(object))[-index]
    ## add a thin argument to the parameter class
    index2 <- index2[seq(1, length(index2), by=thin(param))]
    indices <- sort(c(index, index2))
    object <- object[indices,]
  } else {
    object <- new("MinDistExperiment")
  }
  object
}
