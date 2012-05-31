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
	##trioSet12023chr12 <- updateObjectFromSlots(tSet)
	##genomeBuild(trioSet12023chr12) <- "hg18"
	res <- computeBayesFactor(tSet, gr, prOutlierBAF=list(initial=1e-4, max=1e-2, maxROH=1e-3))
	res <- computeBayesFactor(tSet, rd2, prOutlierBAF=list(initial=1e-4, max=1e-2, maxROH=1e-3))
	checkTrue(as.character(unlist(state(res), use.names=FALSE)) %in% c("334", "333"))
}
