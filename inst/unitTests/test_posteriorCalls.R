test_pennParam <- function(){
  penn <- PennParam()
  checkTrue(validObject(penn))
}

test_MAP2 <- function(){
  library(oligoClasses)
  library(foreach)
  library(BSgenome.Hsapiens.UCSC.hg19)
  foreach::registerDoSEQ()
  data(trioSetListExample)
  me <- as(trioSetList, "MinDistExperiment")
  checkIdentical(me$filename, setNames(c("NA12891", "NA12892", "NA12878"), c("father", "mother", "offspring")))

  seqinfo(me) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(me), ]
  checkTrue(validObject(me))
  me <- subsetAndSort(me, seqlevels(me))
  param <- MinDistParam()
  mdgr <- segment2(me, param)
  mindist(mdgr) <- narrow2(mdgr, param)
  md_g <- unlist(MAP2(me, mdgr, param))
  checkIdentical(sum(md_g$call=="221", na.rm=TRUE),2L)
}

test_posteriorCalls <- function(){
  library(oligoClasses)
  registerDoSEQ()
  gr <- GRanges("chr12", IRanges(21208619, 21520058),
                seg.mean=-0.25,
                mindist.mad=0.19,
                sample="12023_01")
  metadata(gr) <- list(genome="hg18")
  seqlengths(gr) <- getSequenceLengths("hg18")[["chr12"]]
  path <- system.file("extdata", package="MinimumDistance")
  load(file.path(path, "trioSet12023chr12.rda"))
  me <- as(trioSet12023chr12, "MinDistExperiment")
  param <- MinDistParam()
  res <- MAP2(me, gr, param)
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
