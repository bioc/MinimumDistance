setAs("RangedDataCNV", "GRanges", function(from, to){
	cols <- colnames(from)
	gr <- GRanges(paste("chr", chromosome(from), sep=""),
		      IRanges(start(from), end(from)),
		      sample=sampleNames(from),
		      numberProbes=coverage2(from))
	if("seg.mean" %in% cols)
		elementMetadata(gr)$seg.mean <- from$seg.mean
	if("mindist.mad" %in% cols)
		elementMetadata(gr)$mindist.mad <- from$mindist.mad
	gr
})
