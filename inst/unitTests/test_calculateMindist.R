test_calculateMindist <- function(){
	data(trioSetListExample)
	mdlist <- calculateMindist(lrr(trioSetList))

	trioSet <- stack(trioSetList)
	md <- calculateMindist(lrr(trioSet))

	md1 <- do.call("rbind", mdlist)
	dimnames(md1) <- NULL
	dimnames(md) <- NULL
	checkTrue(identical(md, md1))
}
