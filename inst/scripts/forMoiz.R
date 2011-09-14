library(MinimumDistance)
library(oligoClasses)
load("~/Projects/MinimumDistance/data/pedigreeExample-cleft.rda")
load("~/Projects/MinimumDistance/data/samplesheetExample-cleft.rda")
## right now, you must have a column header called 'Sample.Name'
stopifnot("Sample.Name" %in% colnames(samplesheetExample))
samplesheetExample$Sample.Name
## Ideally, the Sample.Name should map back to the filenames of the raw data

##path <- "/thumper/ctsa/beaty/holger/txtfiles"
##filenames <- list.files(path, full.name=TRUE)
cdfName <- "human610quadv1b"

##----------------------------------------------------------------------------
##
## Feature Data
##
##---------------------------------------------------------------------------
fD <- oligoClasses:::featureDataFrom(paste(cdfName, "Crlmm", sep=""))
fD <- fD[order(fD$chromosome, fD$position), ]

##----------------------------------------------------------------------------
##
## Phenotypic Data on trios organized as an array
##
##---------------------------------------------------------------------------
if(FALSE){
	ss <- array(NA, dim=c(nrow(pedigreeExample), ncol(samplesheetExample), 3),
		    dimnames=list(rownames(pedigreeExample),
		    colnames(samplesheetExample),
		    c("F", "M", "O")))
	## in the CleftData, the sample identifier is the first 8 letters of the Sample.Name
	s <- function(x) substr(x, 1, 8)
	father.index <- match(pedigreeExample[, "F"], s(samplesheetExample$Sample.Name))
	mother.index <- match(pedigreeExample[, "M"], s(samplesheetExample$Sample.Name))
	offspring.index <- match(pedigreeExample[, "O"], s(samplesheetExample$Sample.Name))
	ss[, , "F"] <- as.matrix(samplesheetExample[father.index, ])
	ss[, , "M"] <- as.matrix(samplesheetExample[mother.index, ])
	ss[, , "O"] <- as.matrix(samplesheetExample[offspring.index, ])
	stopifnot(all.equal(as.character(s(samplesheetExample$Sample.Name[offspring.index])), as.character(pedigreeExample[, "O"])))
	rownames(ss) <- pedigreeExample[, "O"]
	save(ss, file="~/Projects/MinimumDistance/data/ss-cleft.rda")
} else load("~/Projects/MinimumDistance/data/ss-cleft.rda")


##---------------------------------------------------------------------------
##
##  create an array of BAFs and logR ratios for each chromosome
##
##---------------------------------------------------------------------------
## Hint, to see the markers on chromosome 1.
##chr1.markers <- sampleNames(fD)[fD$chromosome==1]
## Or, create a list of marker names by chromosome
marker.list <- split(sampleNames(fD), fD$chromosome)
##
## I'm going to make a trivial trioSet list with just 100 markers per chromosome
##
marker.list <- lapply(marker.list, sample, size=100)
##
## Let's drop chromosomes 23-25
np <- nrow(pedigreeExample)
marker.list <- marker.list[-(23:25)]
trioSetList <- vector("list", 22)
names(trioSetList) <- 1:22
for(chrom in seq_along(marker.list)){
	## Use the name of the offspring as the name for the trio:
	nr <- length(marker.list[[chrom]])
	bafArray <- logRArray <- array(NA, dim=c(nr, np, 3))
	dimnames(bafArray) <- dimnames(logRArray) <- list(marker.list[[chrom]],
								    as.character(pedigreeExample[, "O"]),
								    c("F", "M", "O"))
	## For each chromosome, create a TrioSet
	pD <- annotatedDataFrameFrom(as.matrix(logRArray[, , 1]), byrow=FALSE)
	sampleNames(pD) <- colnames(logRArray)
	index <- match(marker.list[[chrom]], sampleNames(fD))
	## initialize 'TrioSet'
	trioSetList[[chrom]] <- new("TrioSet",
				    logRRatio=logRArray,
				    BAF=bafArray,
				    phenoData=pD,
				    featureData=fD[index,],
				    mindist=NULL,
				    annotation=cdfName)
	stopifnot(validObject(trioSetList[[chrom]]))
	trioSetList[[chrom]]@phenoData2 <- ss
}
trioSetList <- as(trioSetList, "TrioSetList")
stopifnot(validObject(trioSetList))
save(trioSetList, file="~/Projects/MinimumDistance/data/trioSetListExample-cleft.rda")


example.beadStudioFile <- list.files("/thumper/ctsa/snpmicroarray/hapmap/processed/illumina/610kQuad",
				     pattern="CEU_Final",
				     full.names=TRUE)
c("character", "character", rep(NULL, 7), rep("numeric", 2))
header <- read.csv(example.beadStudioFile, skip=9,
		     colClasses=c("character", "character", rep("NULL", 7), rep("numeric", 2)))
##trace(MinimumDistance:::initializeTrioContainer,browser)
##MinimumDistance::::initializeTrioContainer(pedigree,samplesheet)
##untrace(MinimumDistance:::initializeTrioContainer)


trioRanges <- minimumDistanceCalls(
				   container=tmp,
				   chromosomes=1:22,## for all chromosomes
                                  cbs.filename="ranges.rda") 


md <- minimumDistance(container=trioSetList, readFiles=FALSE, calculate.md=TRUE,
			       container.filename="~/md.rda")

trioSetList <- minimumDistance(container=trioSetList,
                               chromosomes=1:22,
			       readFiles=FALSE,
                               cdfName="human610quadv1b",
                               calculate.md=TRUE,
                               calculate.mad=TRUE)

## how the 'md.segs' object was created
md.segs <- xsegment(object=trioSetList,
                    verbose=TRUE)
## how the 'cbs.segs' object was created
cbs.segs <- xsegment(object=trioSetList,
                     segment.mindist=FALSE,
                     verbose=TRUE)
     
## how the 'called.ranges' object was created
called.ranges <- minimumDistanceCalls(container=trioSetList,
                                      cbs.segs=cbs.segs,
                                      ranges=md.segs,
                                      mindistance.threshold=0.09)



##plot the first third range
             trellisList <- minimumDistancePlot(trioSetList=trioSetList,
                                                ranges=trioRanges[91,],
                                                md.segs=mdRanges,
                                                cbs.segs=cbs.segs
                                                cex=0.2,
                                                scales.cex=0.6)
             logr.trellis.list <- trellisList[[1]]
             baf.trellis.list <- trellisList[[2]]
             gridlayout(lattice.object=list(logr.trellis.list[[1]],
                        baf.trellis.list[[1]]), rd=called.ranges[91,])


