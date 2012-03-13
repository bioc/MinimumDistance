test_cbsSplits <- function(){
	require("human610quadv1bCrlmm")
	path <- system.file("extdata", package="human610quadv1bCrlmm")
	load(file.path(path, "snpProbes.rda"))
	load(file.path(path, "cnProbes.rda"))
	x <- matrix(NA, nrow=nrow(snpProbes)+nrow(cnProbes), 1)
	rownames(x) <- c(rownames(snpProbes), rownames(cnProbes))
	fd <- GenomeAnnotatedDataFrameFrom(x, annotationPkg="human610quadv1bCrlmm")
	fd <- fd[order(chromosome(fd), position(fd)), ]
	fd <- fd[chromosome(fd) < 23, ]
	pos <- position(fd)
	index <- split(seq_len(nrow(fd)), chromosome(fd))
	res <- foreach(i=index) %do% MinimumDistance:::splitByDistance(pos[i], 100e3)
	tabs <- lapply(res, table)
	checkTrue(all(unlist(tabs) > 1000))
}
