\name{xsegment}
\alias{cbs.segs}
\alias{md.segs}
\docType{data}
\title{Segmentation of log R ratios}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
	data(cbs.segs)
	data(md.segs)
}
\format{
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
\dontrun{
        ## show how 'cbs_segs' object was created
	data("trioSetList")
	md.segs <- xsegment(object=trioSetList, verbose=TRUE)
	cbs.segs <- xsegment(object=trioSetList,
			     segment.mindist=FALSE, verbose=TRUE)

	cbs.segs <- vector("list", length(trioSetList))
	chrom <- seq_along(trioSetList)

	for(i in chrom){
		trioSet <- trioSetList[[i]]
		marker.index <- which(!duplicated(position(trioSet)))
		##lR <- logR(trioSet)[marker.index, , ]
		pos <- position(trioSet)[marker.index]
		if(any(diff(order(pos)) < 0)){
			index <- order(pos)
			marker.index <- marker.index[index]
			trioSet <- trioSet[marker.index, ]
			pos <- position(trioSet)
		}
		lR <- logR(trioSet)[marker.index, , , drop=FALSE]
		res <- vector("list", ncol(lR))
		for(j in 1:ncol(lR)){
			CNA.object <- CNA(genomdat=lR[, j, ],
					  chrom=chromosome(trioSet)[marker.index],
					  maploc=pos,
					  data.type="logratio",
					  sampleid=fmoNames(trioSet)[j, ])
			smu.object <- smooth.CNA(CNA.object)
			tmp <- segment(smu.object, verbose=0, alpha=0.01)
			res[[j]] <- cbind(tmp$output, tmp$segRows)
		}
		tmp <- do.call("rbind", res)
		ranges <- IRanges(tmp$loc.start, tmp$loc.end)
		rd.cnv <- RangedDataCBS(ranges=ranges,
					chromosome=tmp$chrom,
					sampleId=tmp$ID,
					coverage=tmp$num.mark,
					seg.mean=tmp$seg.mean,
					startIndexInChromosome=tmp$startRow,
					endIndexInChromosome=tmp$endRow)
		cbs.segs[[i]] <- rd.cnv
	}
	cbs.segs <- stack(RangedDataList(cbs.segs))
	## add the 'trio' name.  Can be the name of the offspring.
	cbs.segs$family <- NA
	father.index <- which(sampleNames(cbs.segs) %in% fatherNames(trioSetList))
	father.names <- sampleNames(cbs.segs)[father.index]
	trio.index <- match(father.names, fatherNames(trioSetList))
	cbs.segs$family[father.index] <- sampleNames(trioSetList)[trio.index]
	mother.index <- which(sampleNames(cbs.segs) %in% motherNames(trioSetList))
	mother.names <- sampleNames(cbs.segs)[mother.index]
	trio.index <- match(mother.names, motherNames(trioSetList))
	cbs.segs$family[mother.index] <- sampleNames(trioSetList)[trio.index]
	offspring.index <- which(sampleNames(cbs.segs) %in% offspringNames(trioSetList))
	offspring.names <- sampleNames(cbs.segs)[offspring.index]
	trio.index <- match(offspring.names, offspringNames(trioSetList))
	cbs.segs$family[offspring.index] <- sampleNames(trioSetList)[trio.index]
	stopifnot(all(trioNames(cbs.segs) %in% sampleNames(trioSetList)))
	save(cbs.segs, file="~/Projects/MinimumDistance/data/cbs.segs.rda")
}
data(cbs.segs)
}
\keyword{datasets}
