# Affy data Minimum Distance Ananlysis
library(MinimumDistance)
library(genomewidesnp6Crlmm)

load("../pedigreeInfo.rda")
load("data.rda")
o.index <- match(colnames(baf), pedigreeInfo[, "O"])
o.index <- o.index[!is.na(o.index)]
#pedigreeExample <- pedigreeExample[1:2, ]
pedigreeInfo <- pedigreeInfo[o.index,]
samplesheetExample <- SampleSheet(filename="NULL",
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
pedigreeExample <- new("Pedigree", trios=pedigreeInfo,
		       trioIndex=pedigreeIndex)
#allNames(pedigreeExample)
#show(pedigreeExample)
trioSetList <- TrioSetList(lrr=logR,
			   baf=baf,
			   pedigree=pedigreeExample,
			   sampleSheet=samplesheetExample,
			   cdfname="genomewidesnp6Crlmm")
md <- calculateMindist(trioSetList)
mindist(trioSetList) <- md
mads.md <- mad2(mindist(trioSetList), byrow=FALSE)
mad.mindist(trioSetList) <- mads.md
mads.lrr.sample <- mad2(lrr(trioSetList), byrow=FALSE)
mads.lrr.marker <- mad2(lrr(trioSetList), byrow=TRUE)
mad.sample(trioSetList) <- mads.lrr.sample
#mad(trioSetList)
mad.marker(trioSetList) <- mads.lrr.marker
#fvarLabels(trioSetList[[1]])
trioSetList <- order(trioSetList)
md.segs <- segment2(object=mindist(trioSetList),
		    pos=position(trioSetList),
		    chrom=chromosome(trioSetList), verbose=1)
#md.segs
#sampleNames(md.segs)
#coverage2(md.segs)
#chromosome(md.segs)
#mean(md.segs)

lrr.segs <- segment2(object=lrr(trioSetList),
		     pos=position(trioSetList),
		     chrom=chromosome(trioSetList),
		     id=trios(trioSetList), ## id is required
		     verbose=TRUE)

md.segs2 <- narrow(md.segs, lrr.segs, 0.9, mad.minimumdistance=mads.md)
nrow(md.segs2) > nrow(md.segs)
map.segs <- computeBayesFactor(object=trioSetList, ranges=md.segs2)
map.segs <- map.segs[order(sampleNames(map.segs), chromosome(map.segs), start(map.segs)), ]
#map.segs
#table(state(map.segs), chromosome(map.segs))

#Visualization
#trioSet <- stack(trioSetList)
#rd.denovoDel <- map.segs[state(map.segs) == 332, ]

#library("VanillaICE")
#lrrfig <- xyplot(lrr ~ x | id,
#		 data=trioSet,
#		 range=rd.denovoDel[1, ],
#		 lrr.segs=lrr.segs,
#		 md.segs=map.segs,
#		 frame=500e3,
#		 panel=xypanelMD,
#                xlab="physical position (Mb)",
#		 cex=0.2,
#		 pch=21,
#		 border="orange",
#                scales=list(cex=0.5),
#		 par.strip.text=list(lines=0.8, cex=0.6),
#		 col.np="royalblue", fill.np="grey",
#		 index.cond=list(4:1),
#		 return.data.frame=FALSE,
#                cex.state=0.5)
#baffig <- xyplot(baf ~ x | id,
#		 data=trioSet,
#		 range=rd.denovoDel[1, ],
#		 xlab="",
#		 lrr.segs=lrr.segs,
#		 md.segs=map.segs,
#		 frame=500e3,
#		 panel=VanillaICE::xypanel,
#		 cex=0.2,
#		 pch=21,
#		 border="orange",
#                scales=list(cex=0.6),
#		 par.strip.text=list(lines=0.8, cex=0.6),
#		 col.np="royalblue", fill.np="grey",
#		 index.cond=list(3:1), cex.state=0.5)

#pdf("MinDistPlot.pdf")
#print(VanillaICE:::arrangeSideBySide(lrrfig, baffig))
#dev.off()

save(trioSetList, lrr.segs, map.segs, file="MinDistData.rda")
rm(list=ls()); gc()

