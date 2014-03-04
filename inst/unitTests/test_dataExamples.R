test_data <- function(){
	data(md.segs)
	if(FALSE){
		md.segs <- oligoClasses:::coerceToGRanges(md.segs, build="hg19")
		save(md.segs, file="~/Software/MinimumDistance/data/md.segs.rda")
	}
	checkTrue(is(md.segs, "GRanges"))

	data(lrr.segs)
	if(FALSE){
		lrr.segs <- oligoClasses:::coerceToGRanges(lrr.segs, build="hg19")
		save(lrr.segs, file="~/Software/MinimumDistance/data/lrr.segs.rda")
	}
	checkTrue(is(lrr.segs, "GRanges"))

	data(map.segs)
	if(FALSE){
		map.segs2 <- oligoClasses:::coerceToGRanges(map.segs, build="hg19")
		elementMetadata(map.segs2)$mindist.mad <- map.segs$mindist.mad
		elementMetadata(map.segs2)$state <- map.segs$state
		elementMetadata(map.segs2)$lik.state <- map.segs$lik.state
		elementMetadata(map.segs2)$lik.norm <- map.segs$lik.norm
		elementMetadata(map.segs2)$argmax <- map.segs$argmax
		map.segs <- map.segs2
		save(map.segs, file="~/Software/MinimumDistance/data/map.segs.rda")
	}
	checkTrue(is(map.segs, "GRanges"))
}
