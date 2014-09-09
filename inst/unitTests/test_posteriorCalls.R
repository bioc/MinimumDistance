test_pennParam <- function(){
  penn <- PennParam()
  checkTrue(validObject(penn))
}

test_MAP2 <- function(){
  library(oligoClasses)
  library(foreach)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(VanillaICE)
  foreach::registerDoSEQ()
  data(trioSetListExample)
  me <- as(trioSetList[,1], "MinDistExperiment")
  seqinfo(me) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(me), ]
  if(FALSE) {
    md_exp <- me
    save(md_exp, file="~/Software/bridge/MinimumDistance/data/md_exp.rda")
  }
  checkTrue(validObject(me))
  ##me <- subsetAndSort(me, seqlevels(me))
  e_param <- EmissionParam(temper=1, p_outlier=1/100)
  penn_param <- PennParam(prNonMendelian=1.5e-6)
  param <- MinDistParam(thin=1L, emission=e_param, penncnv=penn_param)
  mdgr <- segment2(me, param)
  if(FALSE){
    md_gr <- mdgr
    save(md_gr, file="~/Software/bridge/MinimumDistance/data/md_gr.rda")
  }
  md_g <- MAP2(me, mdgr, param)
  checkTrue(length(denovoHemizygous(md_g))==1)

##  findOverlaps(GRanges("chr22", IRanges(20.8e6, 21.4e6)), segs(md_g))
  ##checkIdentical(sum(md_g$call=="221", na.rm=TRUE),2L)
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
  checkTrue(length(denovo(res))==0)
  ##checkTrue(res$call=="222")
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
