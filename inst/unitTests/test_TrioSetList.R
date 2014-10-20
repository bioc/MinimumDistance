test_TrioSetListdataExamples <- function(){
  library(IRanges)
  data(trioSetListExample)
  checkTrue(validObject(trioSetList))
  checkTrue(validObject(trioSetList[[1]]))
  trioSet <- stack(trioSetList)
  checkTrue(validObject(trioSet))
  x <- dim(trioSet)
  y <- setNames(c(43364L, 2L, 3L), c("Features", "Trios", "Members"))
  checkIdentical(x, y)
  library(Biobase)
  checkTrue(storageMode(MinimumDistance:::AssayDataList("lockedEnvironment")) == "lockedEnvironment")
}

test_TrioSet <- function(){
  path <- system.file("extdata", package="MinimumDistance")
  load(file.path(path, "logRratio.rda"))
  load(file.path(path, "baf.rda"))
  load(file.path(path, "pedigreeInfo.rda"))
  ped <- Pedigree(pedigreeInfo)
  trioSet <- TrioSet(lrr=logRratio,
                     baf=baf,
                     pedigree=ped,
                     cdfname="human610quadv1bCrlmm",
                     genome="hg18")
}

test_TrioSetList_construction <- function(){
	library(oligoClasses)
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
				   cdfname="human610quadv1bCrlmm",
				   genome="hg18")
	checkTrue(validObject(trioSetList))
	## the following should fail since hg18 build is not currently
	## available with the supplied annotation package
##	checkException(TrioSetList(lrr=logRratio,
##				   baf=baf,
##				   pedigree=ped,
##				   cdfname="human610quadv1bCrlmm",
##				   genome="hg19"))
	load(file.path(path, "sample.sheet.rda"))
##	checkException(TrioSetList(lrr=logRratio, ## must provide row.names
##				   baf=baf,
##				   pedigree=ped,
##				   sample.sheet=sample.sheet,
##				   cdfname="human610quadv1bCrlmm",
##				   genome="hg18"), silent=TRUE)
	nms <- paste("NA",substr(sample.sheet$Sample.Name, 6, 10),sep="")
	trioSetList <- TrioSetList(lrr=logRratio, ## must provide row.names
				   baf=baf,
				   pedigree=ped,
				   sample.sheet=sample.sheet,
				   row.names=nms,
				   cdfname="human610quadv1bCrlmm",
				   genome="hg18")
	checkTrue(validObject(trioSetList))
	trioSet <- TrioSet(lrr=logRratio,
			   baf=baf,
			   pedigree=ped,
			   cdfname="human610quadv1bCrlmm",
			   genome="hg18")

	if(FALSE){
		path <- system.file("extdata", package="MinimumDistance")
		load(file.path(path, "logRratio.rda"))
		load(file.path(path, "baf.rda"))
		load(file.path(path, "pedigreeInfo.rda"))
		ped <- Pedigree(pedigreeInfo)
		r <- lrr(trioSet)
		chrom <- paste("chr", chromosome(trioSet), sep="")
		seqinfo <- Seqinfo(seqnames=unique(chrom),
				   genome="hg18")
		rowData <- GRanges(chrom,
				   IRanges(position(trioSet)-12,
					   position(trioSet)+12),
				   seqinfo=seqinfo)
		## rownames of colData are the column names of the
		## summarized experiment - rownames are offspring ids.
		## - colnames are father, mother, offspring - 3rd
		## dimension is covariates on the family members
		## (e.g., age, ...)  - perhaps require id for a label
		## in third dimension.
		## class Pedigree could inherit from DataFrame.
		ary <- I(array(NA, dim=c(ncol(trioSet), 3, 2)))
		dimnames(ary) <- list(sampleNames(ped),
				      c("father", "mother", "offspring"),
				      c("id", "age"))
		##colData <- DataFrame(ary, metadata=rep(c("father", "mother", "offspring"), 2))
		colData <- DataFrame(ary) ##metadata=rep(c("father", "mother", "offspring"), 2))
		## Need to rewrite fatherNames, motherNames, offspringNames accessors
		colData[[1]][, "father", "id"] ##make an accessor
		colData[[1]][, "mother", "id"] ##make an accessor
		colData[[1]][, "offspring", "id"] ##make an accessor
		rownames(colData) <- sampleNames(ped)
		##colnames(colData) <- c("father", "mother", "offspring")
		se <- SummarizedExperiment(assays=SimpleList(lrr=r),
					   rowData=rowData,
					   colData=colData)
	}
	checkTrue(validObject(trioSet))
	trioSet <- TrioSet(lrr=logRratio,
			   baf=baf,
			   pedigree=ped,
			   sample.sheet=sample.sheet,
			   row.names=nms,
			   cdfname="human610quadv1bCrlmm",
			   genome="hg18")
	checkTrue(validObject(trioSet))

	checkTrue(validObject(trioSetList[[1]]))
	obj <- trioSetList[c(1,2)]
	checkIdentical(chromosome(obj), c(1L, 2L))
	obj <- trioSetList[c(1,1)]
	checkIdentical(chromosome(obj), c(1L, 1L))
	obj <- trioSetList[FALSE]
	checkTrue(validObject(obj))
	checkIdentical(length(obj), 0L)

	b <- baf(trioSetList)
	b <- b[-1]
	library(Biobase)
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
	checkTrue(is(trioSet, "TrioSet"))

	checkTrue(validObject(trioSet[1, ]))
	triosubset <- trioSet[1:5, 1]
	checkIdentical(as.integer(dim(triosubset)), c(5L, 1L, 3L))
	triosubset <- trioSet[1:5, ]
	checkIdentical(as.integer(dim(triosubset)), c(5L, 2L, 3L))
	triosubset <- trioSet[, 1]
	checkIdentical(as.integer(dim(triosubset)), c(as.integer(nrow(trioSet)), 1L, 3L))
}


##
##test_TrioSetListLD <- function(){
##  ## constructor for large data
##  library(oligoClasses)
##  path <- system.file("extdata", package="MinimumDistance")
##  fnames <- list.files(path, pattern=".txt")
##  ##allow duplicated father and mother names
##  ped <- Pedigree(data.frame(F=c("F.txt", "F.txt"),
##                             M=c("M.txt", "M.txt"),
##                             O=c("O.txt", "O1.txt")))
##
##  ##datlist <- lapply(fnames, VanillaICE::read.bsfiles)
##
##  ##dat <- VanillaICE::read.bsfiles(file.path(path, fnames))
##  ##trace(TrioSetListLD, browser)
##  trioSetList <- TrioSetListLD(path=path,
##                               fnames=fnames,
##                               pedigreeData=ped,
##                               annotationPkg="human610quadv1bCrlmm",
##                               genome="hg18")
##  checkTrue(validObject(trioSetList))
##  checkTrue(is(lrr(trioSetList)[[1]], "array"))
##
##  library2(ff)
##  library2(foreach)
##  ldPath(tempdir())
##  registerDoSEQ()
##  trioSetListff <- TrioSetListLD(path=path,
##                                 fnames=fnames,
##                                 pedigreeData=ped,
##                                 annotationPkg="human610quadv1bCrlmm",
##                                 genome="hg18")
##  checkTrue(validObject(trioSetListff))
##  checkTrue(is(lrr(trioSetListff)[[1]], "ff_array"))
##  checkTrue(identical(lrr(trioSetListff)[[1]][,,], lrr(trioSetList)[[1]]))
##}
