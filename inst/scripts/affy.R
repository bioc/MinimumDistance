#Affy HapMap data
library(ff)
ldPath("/amber1/archive/oncbb/rscharpf/MinDistCompendium/inst/intdata")
library(crlmm)
load("/amber1/archive/oncbb/rscharpf/MinDistCompendium/data/cnSet_hapmapIII.rda")
tmpdir <- "/amber1/archive/oncbb/rscharpf/MinDistCompendium/tmpFiles"
#keep only autosomes
#dim(cnSet)
snames <- sampleNames(cnSet)
# Pedigree information for all HapMap samples
rels <- read.table("/home/mmi/mbootwal/projects/packages/PennCNV/relationships_w_pops_121708.txt",
		   header=TRUE, stringsAsFactors=FALSE)
snames <- snames[snames%in%rels$IID]
# find trios from rels
trio <- cbind(as.integer(rels$IID%in%snames), as.integer(rels$dad%in%snames), as.integer(rels$mom%in%snames))
sum <- apply(trio, 1, sum)
rels2 <- rels[sum == 3, ]
rels2.names <- unlist(rels2[, 2:4])
# index of samples in trios
keep.index <- which(sampleNames(cnSet) %in% rels2.names)
# assign batch name to each sample
batch <- cbind(sampleNames(cnSet)[keep.index], batch(cnSet)[keep.index])
# create a family matrix
fam <- cbind("F"=rels2$dad, "M"=rels2$mom, "O"=rels2$IID)
# Add batch information for each family member
f.index <- match(fam[,1], batch[,1])
m.index <- match(fam[,2], batch[,1])
o.index <- match(fam[,3], batch[,1])
fam <- cbind(fam, "fBatch"=batch[f.index,2], "mBatch"=batch[m.index,2], "oBatch"=batch[o.index,2])
# Remove trios from batches CHEAP, CORER, TESLA
i <- c(which(fam[,4] == "CORER" | fam[,5] == "CORER" | fam[,6] == "CORER"), which(fam[,4] == "CHEAP" | fam[,5] == "CHEAP" | fam[,6] == "CHEAP"), which(fam[,4] == "TESLA" | fam[,5] == "TESLA" | fam[,6] == "TESLA"))
fam <- fam[-i,]
# Index each trio member in order to sampleNames of cnSet object in order to subset it
##tNames <- list()
##for(i in seq_along(fam[,1])){
##	tmp <- unlist(as.character(fam[i,1:3]))
##	tNames[[i]] <- tmp
##}
##tNames <- unlist(tNames)
tNames <- as.character(fam[, 1:3])
s.index <- match(tNames, sampleNames(cnSet))
batch.names <- unique(batch(cnSet)[s.index])
f.index <- split(seq_len(nrow(cnSet)), chromosome(cnSet))
if(FALSE){
	for(i in seq_len(22)){
		message("chromosome ", i)
		logR.baf <- list()
		tmpSet <- cnSet[f.index[[i]], s.index]
		for(j in seq_along(batch.names)){
			B <- batch.names[j]
			logR.baf[[j]] <- calculateRBaf(tmpSet, B)
		}
		lrrs <- lapply(logR.baf, "[[", i="lrr")
		lrrs <- do.call("cbind", lrrs)
		bafs <- lapply(logR.baf, "[[", i="baf")
		bafs <- do.call("cbind", bafs)
		save(lrrs, file=file.path(tmpdir, paste("lrrs_chrom", i, ".rda", sep="")))
		save(bafs, file=file.path(tmpdir, paste("bafs_chrom", i, ".rda", sep="")))
		rm(tmpSet, bafs, lrrs, logR.baf); gc()
	}
}
q('no')
##chrom <- which(chromosome(cnSet) <= 22 & !is.na(chromosome(cnSet)))
##cnSet.trios <- cnSet[chrom, index]
lrr.fns <- list.files(tmpdir, pattern="lrrs_chrom", full.names=TRUE)
baf.fns <- list.files(tmpdir, pattern="bafs_chrom", full.names=TRUE)

load(baf.fns[1])
load("/amber1/archive/oncbb/rscharpf/moiz/Affy/pedigreeInfo.rda")
##load("data.rda")
o.index <- match(colnames(bafs), pedigreeInfo[, "O"])
o.index <- o.index[!is.na(o.index)]
#pedigreeExample <- pedigreeExample[1:2, ]
pedigreeInfo <- pedigreeInfo[o.index,]
samplesheetHapmap <- SampleSheet(filename="NULL",
				 id=as.character(unlist(pedigreeInfo)),
				 plate="NULL")
#pedigreeInfo <- pedigreeExample
#rownames(pedigreeInfo) <- NULL
indId <- unlist(pedigreeInfo)
trioFactor <- matrix(seq_len(nrow(pedigreeInfo)), nrow(pedigreeInfo), 3, byrow=FALSE)
memberId <- substr(names(indId), 1, 1)
## which row in the pedigree file
trio.index <- as.integer(trioFactor) ## note this gives the index in the pedigree
pedigreeIndex <- data.frame(individualId=indId,
			    memberId=memberId,
			    index.in.pedigree=trio.index,
			    stringsAsFactors=FALSE)
rownames(pedigreeIndex) <- NULL
pedigreeHapmap <- new("Pedigree", trios=pedigreeInfo,
		       trioIndex=pedigreeIndex)

## read in the files and create a TrioSetList object for each
## chromosome (length 1), then combine the TrioSetList objects into a
## single TrioSetList
library(MinimumDistance)
##for(i in seq_along(lrr.fns)){
res <- vector("list", length(baf.fns))

for(i in 2:length(lrr.fns)){
	load(lrr.fns[i])
	load(baf.fns[i])
	chr <- strsplit(lrr.fns[i], "chrom")[[1]][[2]]
	chr <- strsplit(chr, ".rda")[[1]][[1]]
	res[[i]] <- TrioSetList(lrr=lrrs,
				baf=bafs,
				pedigree=pedigreeHapmap,
				sampleSheet=samplesheetHapmap,
				cdfname="genomewidesnp6Crlmm",
				chromosome=chr)
	names(res)[i] <- chr
}
## create a single trioSetList
trioSetList <- new("TrioSetList",
		   pedigreeData=pedigreeHapmap,
		   sampleSheet=samplesheetHapmap)
for(i in seq_along(res)) trioSetList[[i]] <- res[[i]][[1]]
names(trioSetList) <- names(res)
index <- order(as.numeric(names(res)))
trioSetList <- trioSetList[index]
names(trioSetList) <- names(res)[index]
save(trioSetList, file="/amber1/archive/oncbb/rscharpf/MinDistCompendium/data/trioSetList.rda")

if(FALSE){
	load("/amber1/archive/oncbb/rscharpf/MinDistCompendium/data/trioSetList.rda")
	## test on 2 chromosomes and 2 trios
	trioSetList <- trioSetList[10:11]
	trace("[", browser, signature="TrioSetList")
	trioSetList <- trioSetList[, 1:2]
	save(trioSetList, file="/amber1/archive/oncbb/rscharpf/MinDistCompendium/tmpFiles/trioSetList_small.rda")
}

library(MinimumDistance)
library(ff)
load("/amber1/archive/oncbb/rscharpf/MinDistCompendium/tmpFiles/trioSetList_small.rda")
md <- calculateMindist(trioSetList)
mindist(trioSetList) <- md
mads.md <- mad2(mindist(trioSetList), byrow=FALSE)
mad.mindist(trioSetList) <- mads.md
mads.lrr.sample <- mad2(lrr(trioSetList), byrow=FALSE)
mads.lrr.marker <- tryCatch(mad2(lrr(trioSetList), byrow=TRUE), error=function(e) NULL)
if(is.null(mads.lrr.marker)){
	mads.lrr.marker <- list()
	mads.lrr.marker[[1]] <- rep(mads.lrr.sample[1, "O"], nrow(trioSetList[[1]]))
	mads.lrr.marker[[2]] <- rep(mads.lrr.sample[2, "O"], nrow(trioSetList[[2]]))
}
mad.sample(trioSetList) <- mads.lrr.sample
mad(trioSetList)
mad.marker(trioSetList) <- mads.lrr.marker
fvarLabels(trioSetList[[1]])
md.segs <- segment2(object=mindist(trioSetList),
		    pos=position(trioSetList),
		    chrom=chromosome(trioSetList), verbose=1)

lrr.segs <- segment2(object=lrr(trioSetList),
		     pos=position(trioSetList),
		     chrom=chromosome(trioSetList),
		     id=trios(trioSetList), ## id is required
		     verbose=TRUE)

md.segs2 <- narrow(md.segs, lrr.segs, 0.9, mad.minimumdistance=mads.md)
map.segs <- computeBayesFactor(object=trioSetList, ranges=md.segs2)
map.segs <- map.segs[order(sampleNames(map.segs), chromosome(map.segs), start(map.segs)), ]

trioSet <- stack(trioSetList)
##rd.denovoDel <- map.segs[state(map.segs) == 332, ]
library("VanillaICE")
lrrfig <- xyplot(lrr ~ x | id,
		 data=trioSet,
		 range=map.segs[7, ],
		 lrr.segs=lrr.segs,
		 md.segs=map.segs,
		 frame=500e3,
		 panel=xypanelMD,
                 xlab="physical position (Mb)",
		 cex=0.2,
		 pch=21,
		 border="orange",
                 scales=list(cex=0.5),
		 par.strip.text=list(lines=0.8, cex=0.6),
		 col.np="royalblue", fill.np="grey",
		 index.cond=list(4:1),
		 return.data.frame=FALSE,
                 cex.state=0.5)
baffig <- xyplot(baf ~ x | id,
		 data=trioSet,
		 range=map.segs[7, ],
		 xlab="",
		 frame=500e3,
		 panel=VanillaICE::xypanel,
		 cex=0.2,
		 pch=21,
		 border="orange",
                 scales=list(cex=0.6),
		 par.strip.text=list(lines=0.8, cex=0.6),
		 col.np="royalblue", fill.np="grey",
		 index.cond=list(3:1), cex.state=0.5)
library(VanillaICE)
print(arrangeSideBySide(lrrfig, baffig))
