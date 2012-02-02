test_dataExamples <- function(){
	data(trioSetListExample)
	checkTrue(validObject(trioSetList))
	trioSet <- stack(trioSetList)
	checkTrue(validObject(trioSet))
}

test_TrioSetList_construction <- function(){
	checkTrue(validObject(new("TrioSetList")))
	checkTrue(validObject(TrioSetList()))
	checkTrue(validObject(TrioSetList(chromosome=1:22)))
	checkException(TrioSetList(chromosome=1:23), silent=TRUE)

	path <- system.file("extdata", package="MinimumDistance")
	load(file.path(path, "logRratio.rda"))
	load(file.path(path, "baf.rda"))
	load(file.path(path, "pedigreeInfo.rda"))
	ped <- Pedigree(pedigreeInfo)
	trioSetList <- TrioSetList(lrr=logRratio,
				   baf=baf,
				   pedigree=ped,
				   cdfname="human610quadv1bCrlmm")
	checkTrue(validObject(trioSetList))
	load(file.path(path, "sample.sheet.rda"))
	checkException(TrioSetList(lrr=logRratio, ## must provide row.names
				   baf=baf,
				   pedigree=ped,
				   sample.sheet=sample.sheet,
				   cdfname="human610quadv1bCrlmm"), silent=TRUE)
	nms <- paste("NA",substr(sample.sheet$Sample.Name, 6, 10),sep="")
	trioSetList <- TrioSetList(lrr=logRratio, ## must provide row.names
					  baf=baf,
					  pedigree=ped,
					  sample.sheet=sample.sheet,
					  row.names=nms,
					  cdfname="human610quadv1bCrlmm")
	checkTrue(validObject(trioSetList))
	trioSet <- TrioSet(lrr=logRratio,
			   baf=baf,
			   pedigree=ped,
			   cdfname="human610quadv1bCrlmm")
	checkTrue(validObject(trioSet))
	trioSet <- TrioSet(lrr=logRratio,
			   baf=baf,
			   pedigree=ped,
			   sample.sheet=sample.sheet,
			   row.names=nms,
			   cdfname="human610quadv1bCrlmm")
	checkTrue(validObject(trioSet))

	checkTrue(validObject(trioSetList[[1]]))
	obj <- trioSetList[c(1,2)]
	checkIdentical(chromosome(obj), c(1L, 2L))
	obj <- trioSetList[c(1,1)]
	checkIdentical(chromosome(obj), c(1L, 1L))
	obj <- trioSetList[FALSE]
	checkTrue(validObject(obj))

	b <- baf(trioSetList)
	b <- b[-1]
	object <- assayDataElementReplace(trioSetList, "BAF", b)
	checkException(validObject(object), silent=TRUE)
	b <- baf(trioSetList)
	b[[1]] <- b[[1]][, , 1:2]
	object <- assayDataElementReplace(trioSetList, "BAF", b)
	checkException(validObject(object), silent=TRUE)

	object <- trioSetList
	object@chromosome <- chromosome(trioSetList)[1]
	checkException(validObject(object), silent=TRUE)

	object <- trioSetList
	object@featureDataList <- object@featureDataList[1:2]
	checkException(validObject(object), silent=TRUE)

	## TrioSet construction
	checkTrue(validObject(new("TrioSet")))
	trioSet <- trioSetList[[1]]
	checkTrue(validObject(trioSet))
	checkTrue(is(trioSet, "TrioSet"))

	checkTrue(validObject(trioSet[1, ]))
	triosubset <- trioSet[1:5, 1]
	checkIdentical(as.integer(dim(triosubset)), c(5L, 1L, 3L))
	triosubset <- trioSet[1:5, ]
	checkIdentical(as.integer(dim(triosubset)), c(5L, 2L, 3L))
	triosubset <- trioSet[, 1]
	checkIdentical(as.integer(dim(triosubset)), c(25L, 1L, 3L))

}



test_TrioSetListLD <- function(){
	## constructor for large data
	path <- system.file("extdata", package="MinimumDistance")
	fnames <- list.files(path, pattern=".txt")
	##allow duplicated father and mother names
	ped <- Pedigree(data.frame(F=c("F.txt", "F.txt"),
				   M=c("M.txt", "M.txt"),
				   O=c("O.txt", "O1.txt")))
	trioSetList <- TrioSetListLD(path=path,
				     fnames=fnames,
				     pedigreeData=ped,
				     annotationPkg="human610quadv1bCrlmm")
	checkTrue(validObject(trioSetList))
	checkTrue(is(lrr(trioSetList)[[1]], "array"))

	library(ff)
	ldPath(tempdir())
	trioSetListff <- TrioSetListLD(path=path,
				       fnames=fnames,
				       pedigreeData=ped,
				       annotationPkg="human610quadv1bCrlmm")
	checkTrue(validObject(trioSetListff))
	checkTrue(is(lrr(trioSetListff)[[1]], "ff_array"))
	checkTrue(identical(lrr(trioSetListff)[[1]][,,], lrr(trioSetList)[[1]]))
}


test_subsetTrioSetList <- function(){
	object <- Pedigree(fatherIds=rep(letters[1:3], each=2),
			   motherIds=rep(letters[4:6], each=2),
			   offspringIds=letters[11:16])
	checkTrue(validObject(object[6, ]))
}
