make_test_Pedigree <- function(){
	pedTest <- new("Pedigree", trios=data.frame(F=c("NA06993", "NA11881"),
				   M=c("NA06985", "NA11882"),
				   O=c("NA06991", "NA10859"),
				   stringsAsFactors=FALSE),
		       trioIndex=data.frame(individualId=c("NA06993", "NA11881",
				   "NA06985", "NA11882", "NA06991", "NA10859"),
		       memberId=rep(c("F", "M", "O"), each=2),
		       index.in.pedigree=rep(1:2, 3),
		       stringsAsFactors=FALSE))
}

test_Pedigree_construction <- function(){
	##checkException(Pedigree(), silent=TRUE)
	##trace(Pedigree, browser)
	checkTrue(validObject(Pedigree()))
	path <- system.file("extdata", package="MinimumDistance")
	load(file.path(path, "pedigreeInfo.rda"))
	ped <- Pedigree(fatherIds=pedigreeInfo$F,
			motherIds=pedigreeInfo$M,
			offspringIds=pedigreeInfo$O)
	checkTrue(validObject(ped))
	checkTrue(validObject(Pedigree(pedigreeInfo)))
	ped2 <- make_test_Pedigree()
	checkIdentical(ped, ped2)

	## can not have duplicate offspring identifiers
	checkException(Pedigree(data.frame(F=c("F0.txt", "F0.txt"),
					   M=c("M0.txt", "M0.txt"),
					   O=c("O0.txt", "O0.txt"))), silent=TRUE)

	checkTrue(validObject(Pedigree(data.frame(F=c("F0.txt", "F0.txt"),
						  M=c("M0.txt", "M0.txt"),
						  O=c("O0.txt", "O1.txt")))))

	trio <- trios(ped)
	trio$F[[1]] <- "badname"
	ped@trios <- trio
	checkException(validObject(ped), silent=TRUE)
}
