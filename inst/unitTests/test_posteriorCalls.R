test_pipeline <- function(){
	library(MinimumDistance); library(RUnit)
	data(trioSetListExample)
	md <- calculateMindist(lrr(trioSetList))
	md.segs <- segment2(trioSetList, md=md)
	metadata(md.segs)
	lrr.segs <- segment2(trioSetList, segmentParents=TRUE)
	metadata(lrr.segs)
	mads.md <- mad2(md, byrow=FALSE)
	md.segs2 <- narrow(md.segs, lrr.segs, thr=0.75, mad.minimumdistance=mads.md, fD=featureData(trioSetList))
	map.segs <- computeBayesFactor(object=trioSetList, ranges=md.segs2, prOutlierBAF=list(initial=1e-4, max=1e-2, maxROH=1e-3))
}

test_posteriorCalls <- function(){
	library(MinimumDistance); library(RUnit)
	library(oligoClasses)
	library2(foreach)
	registerDoSEQ()
	library(GenomicRanges)
	library2(MinimumDistance)
	library2(Biobase)
	gr <- GRanges("chr12", IRanges(21208619, 21520058),
		      seg.mean=-0.25,
		      mindist.mad=0.19,
		      sample="12023_01")
	seqlengths(gr) <- getSequenceLengths("hg18")[["chr12"]]
	rd2 <- RangedDataCNV(IRanges(21208619, 21520058),
			    sampleId="12023_01",
			    chrom=12,
			    seg.mean=-0.25,
			    mindist.mad=0.19)
	path <- system.file("extdata", package="MinimumDistance")
	load(file.path(path, "trioSet12023chr12.rda"))
	tSet <- trioSet12023chr12
	res <- computeBayesFactor(tSet, gr, prOutlierBAF=list(initial=1e-4, max=1e-2, maxROH=1e-3))
	res <- computeBayesFactor(tSet, rd2, prOutlierBAF=list(initial=1e-4, max=1e-2, maxROH=1e-3))
	##checkTrue(as.character(unlist(state(res), use.names=FALSE)) %in% c("334", "333"))
	if(FALSE){
		ylab <- expression(log[2]("R ratios"))
		ylab2 <- expression("B allele freqencies")
		at <- c(-1, 0, log2(3/2), log2(4/2))
		labels <- expression(-1, 0, log[2](3/2), log[2](4/2))
		fig <- xyplotTrio(rd=res,
				  object=tSet,
				  frame=2e6,
				  ylab="",
				  xlab="physical position (Mb)",
				  panel=xypanelTrio,
				  scales=list(cex=0.7, x=list(relation="same"),
				  y=list(alternating=1, at=at, labels=labels)),
				  ##lrr.segments=rsegs,
				  ##md.segments=msegs,
				  col.hom="grey50",
				  col.het="grey50",
				  col.np="grey10",
				  segment.col="black",
				  state.cex=0.8,
				  pch=".",
				  cex=1,
				  layout=c(1, 4),
				  ylim=c(-4, 2),
				  key=list(text=list(c(ylab, ylab2),
					   col=c("grey50", "blue")), columns=2))
	}
}
