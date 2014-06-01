test_pennParam <- function(){
  penn <- PennParam()
  checkTrue(validObject(penn))
}

test_MAP2 <- function(){
  library(oligoClasses)
  library(foreach)
  foreach::registerDoSEQ()
  data(trioSetListExample)
  md <- calculateMindist(lrr(trioSetList))
  path <- tryCatch(system.file("unitTests", package="MinimumDistance", mustWork=TRUE),
                   error=function(e) "~/Software/bridge/MinimumDistance/inst/unitTests")
  if(FALSE){
    md.segs <- segment2(trioSetList, md=md, verbose=0)
    saveRDS(md.segs, file=file.path(path, "md_segs.rds"))
  } else md.segs <- readRDS(file.path(path, "md_segs.rds"))
  if(FALSE){
    lrr.segs <- segment2(trioSetList, segmentParents=TRUE, verbose=0)
    saveRDS(lrr.segs, file=file.path(path,"lrr_segs.rds"))
  } else lrr.segs <- readRDS(file.path(path, "lrr_segs.rds"))
  ##metadata(lrr.segs)
  mads.md <- mad2(md, byrow=FALSE)
  ##trace(narrowRanges, browser)
  md.segs2 <- sort(narrowRanges(md.segs, lrr.segs, thr=0.75,
                                mad.minimumdistance=mads.md,
                                fD=featureData(trioSetList)))
  trioSet <- stack(trioSetList)
  se <- as(trioSet, "SnpArrayExperiment")
  rownames(se) <- featureNames(trioSet)
  genome(se) <- "hg19"
  ##  library(devtools)
  ##load_all("~/Software/bridge/VanillaICE")
  ##trace(computeEmissionProbs, browser)
##  E <- computeEmissionProbs(se)
  param <- PennParam()
##  checkException(compute_loglik(E, md_ranges=md.segs2, param=penn))
##  md_ranges <- sort(md.segs2[md.segs2$sample %in% colnames(se)[3]])
##  untrace(compute_loglik, browser)
##  md_ranges <- compute_loglik(E, md_ranges=md_ranges, param=param)
  md_ranges <- MAP2(se, md.segs2, param)
  checkTrue(sum(md_ranges$call=="221", na.rm=TRUE)==3)
##
##  ##
##  ## 1. Make MAP work with new VI interface by coercing TrioSet to a SnpArrayExperiment
##  ## 2. Repeat for TrioSetList
##  ## 3. Replace TrioSets with TrioExperiment classes
##  ## 4. Deprecate old classes
##  ## map.segs <- MAP(trioSet, md.segs2)
}

test_posteriorCalls <- function(){
  library(oligoClasses)
  ##library2(foreach)
  registerDoSEQ()
  ##library2(Biobase)
  gr <- GRanges("chr12", IRanges(21208619, 21520058),
                seg.mean=-0.25,
                mindist.mad=0.19,
                sample="12023_01")
  metadata(gr) <- list(genome="hg18")
  seqlengths(gr) <- getSequenceLengths("hg18")[["chr12"]]
  path <- system.file("extdata", package="MinimumDistance")
  load(file.path(path, "trioSet12023chr12.rda"))
  tSet <- trioSet12023chr12
  ##res <- MAP(tSet, gr, prOutlierBAF=list(initial=1e-4, max=1e-2, maxROH=1e-3))
  se <- as(tSet, "SnpArrayExperiment")
  param <- PennParam(minimum_distance_threshold=0.1)
  gr$sample <- colnames(se)[3]
  res <- MAP2(se, gr, param)
  checkTrue(res$call=="222")
  ##checkTrue(as.character(unlist(state(res), use.names=FALSE)) %in% c("334", "333"))
##  if(FALSE){
##    ylab <- expression(log[2]("R ratios"))
##    ylab2 <- expression("B allele freqencies")
##    at <- c(-1, 0, log2(3/2), log2(4/2))
##    labels <- expression(-1, 0, log[2](3/2), log[2](4/2))
##    fig <- xyplotTrio(rd=res,
##                      object=tSet,
##                      frame=2e6,
##                      ylab="",
##                      xlab="physical position (Mb)",
##                      panel=xypanelTrio,
##                      scales=list(cex=0.7, x=list(relation="same"),
##                        y=list(alternating=1, at=at, labels=labels)),
##                      col.hom="grey50",
##                      col.het="grey50",
##                      col.np="grey10",
##                      segment.col="black",
##                      state.cex=0.8,
##                      pch=".",
##                      cex=1,
##                      layout=c(1, 4),
##                      ylim=c(-4, 2),
##                      key=list(text=list(c(ylab, ylab2),
##                                 col=c("grey50", "blue")), columns=2))
##  }
}
