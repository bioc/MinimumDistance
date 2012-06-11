xyplotTrio <- function(rd, object, frame=200e3, lrr.segments=NULL, md.segments=NULL, nchar_sampleid=15L, ...){
	##if(!is(rd, "RangedDataCNV")) stop("rd is not a RangedDataCNV-derived class")
	if(!is(object, "TrioSet")) stop("object is not a TrioSet")
	if(is.null(mindist(object))) stop("must add minimum distance matrix to mindist slot. Use mindist(object) <- value")
	index <- seq_len(length(rd))
	i <- NULL
	df <- foreach(i=index, .combine="rbind") %do% {
		dataFrameFromRange2(range=rd[i, ],
				    object=object,
				    frame=frame,
				    range.index=i,
				    nchar_sampleid=nchar_sampleid)
	}
	df$range <- factor(paste("range", df$range), ordered=TRUE, levels=unique(paste("range", df$range)))
	index <- split(seq_len(nrow(df)), df$range)
	i <- j <- NULL
	figs <- foreach(i=index, j=seq_along(index)) %do% {
		xyplot(y~x|memberId,
		       data=df[i, ],
		       baf=df$baf[i],
		       is.snp=df$is.snp[i],
		       range=rd[j, ],
		       memberId=df$memberId[i],
		       lrr.segments=lrr.segments,
		       md.segments=md.segments,
		       ped=pedigree(object), ...)
	}
	if(length(figs)==1) figs <- figs[[1]]
	return(figs)
}


##xyplotTrioListLrrBaf <- function(rd, md, object, frame=200e3,
##				 lrr.segments, md.segments, ...){
##	## assume rd is one range
##	object <- object[[chromosome(rd)]]
##	marker.index <- subjectHits(findOverlaps(rd, featureData(object), maxgap=frame))
##	trio.index <- match(sampleNames(rd), sampleNames(object))
##	object <- object[marker.index, trio.index]
##	md <- md[[chromosome(rd)]]
##	md <- md[marker.index, trio.index, drop=FALSE]
##	mindist(object) <- md
## 	xyplotTrioLrrBaf(rd=rd, object=object, frame=frame, lrr.segments=lrr.segments, md.segments=md.segments, ...)
##}


##xyplotTrioSetList <- function(object,
##			      mdlist,
##			      frame=200e3,
##			      map.segment,
##			      lrr.segments,
##			      md.segments, ...){
##	map.segment <- map.segment[1,]
##	CHR <- chromosome(map.segment)
##	id <- ss(sampleNames(map.segment))
##	lrr.segments <- lrr.segments[chromosome(lrr.segments) == CHR & ss(sampleNames(lrr.segments)) == id, ]
##	md.segments <- md.segments[chromosome(md.segments) == CHR & ss(sampleNames(md.segments)) == id, ]
##	ix <- match(sampleNames(map.segment), sampleNames(object))
##	trioSet <- object[[chromosome(map.segment)]][, ix]
##	mindist(trioSet) <- mdlist[[chromosome(map.segment)]][,ix]
##	ylab <- expression(log[2]("R ratios"))
##	ylab2 <- expression("B allele frequencies")
##	at <- c(-1, 0, log2(3/2), log2(4/2))
##	labels <- expression(-1, 0, log[2](3/2), log[2](4/2))
##	df <- dataFrameFromRange2(range=map.segment,
##				  object=trioSet,
##				  frame=frame,
##				  range.index=1)
##	df$range <- factor(paste("range", df$range), ordered=TRUE, levels=unique(paste("range", df$range)))
##	fig <- xyplot(y~x|memberId,
##		      data=df,
##		      baf=df$baf,
##		      is.snp=df$is.snp,
##		      range=map.segment,
##		      memberId=df$memberId,
##		      lrr.segments=lrr.segments,
##		      md.segments=md.segments,
##		      ped=pedigree(object),
##		      ylab="",
##		      xlab="physical position (Mb)",
##		      panel=xypanelTrio,
##		      scales=list(x="same",
##		      y=list(alternating=1, at=at, labels=labels)),
##		      layout=c(1, 4),
##		      ylim=c(-3, 1.5),
##		      key=list(text=list(c(ylab, ylab2),
##			       col=c("black", "blue")), columns=2),
##		      ...)
##	results <- list(trellis=fig,
##			trioSet=trioSet,
##			map.segment=map.segment,
##			lrr.segments=lrr.segments,
##			md.segments=md.segments)
##	return(results)
##}

xypanelTrio <- function(x, y,
			memberId,
			baf,
			is.snp,
			range,
			lrr.segments,
			md.segments,
			baf.color="blue",
			col.hom="grey20",
			fill.hom="lightblue",
			col.het="grey20" ,
			fill.het="salmon",
			col.np="grey20",
			fill.np="grey60",
			state.show=TRUE,
			state.cex=1,
			state.col="blue",
			segment.col="grey50",
			ped,
			..., subscripts){
	panel.abline(h=c(-1, 0, log2(3/2), log2(4/2)), col="grey", lty=2)
	panel.xyplot(x[1], y[1], col="white", ...) ## set it up, but don't plot
	is.snp <- is.snp[subscripts]
	ylim <- current.panel.limits()$ylim
	y[y>ylim[2]] <- ylim[2]
	##
	lpoints(x[!is.snp], y[!is.snp], col=col.np, fill=fill.np, ...)
	## use whatever col.hom to color SNPs
	lpoints(x[is.snp], y[is.snp], col=col.hom, fill=fill.hom, ...)
	j <- panel.number()
	st <- start(range)[j]/1e6
	lrect(xleft=st, xright=end(range)[j]/1e6, ybottom=-10, ytop=10, ...)
	if(state.show){
		## left justify the label to the start of the range
		y.max <- ylim[2]
		ltext(st, y.max, labels=paste("state", state(range)[j]),
		      adj=c(0,1), cex=state.cex, col=state.col)
	}
	b <- baf[subscripts]
	b[b==2] <- NA
	blim <- c(ylim[1]+0.1, ylim[1]+1.5)
	bnew <- rescale(b, blim[1], blim[2])
	lpoints(x[is.snp], bnew[is.snp], col=baf.color, ...)
	memberId <- unique(memberId[subscripts])
	mindistPanel <- match(memberId, levels(memberId))==1
	##sns <- unique(sampleNames[subscripts])
	if(mindistPanel){
		if(!is.null(md.segments)){
			md.segments <- md.segments[sampleNames(md.segments) %in% sampleNames(range) & as.character(chromosome(md.segments)) == as.character(chromosome(range)), ]
			lsegments(x0=start(md.segments)/1e6,
				  x1=end(md.segments)/1e6,
				  y0=elementMetadata(md.segments)$seg.mean,
				  y1=elementMetadata(md.segments)$seg.mean, lwd=2, col=segment.col)		}
	} else {
		## range is labeled by offspring id.
		if(!is.null(lrr.segments)){
			j <- match(sampleNames(range), sampleNames(ped))
			fatherPanel <- match(memberId, levels(memberId))==4
			motherPanel <- match(memberId, levels(memberId))==3
			offsprPanel <- match(memberId, levels(memberId))==2
			if(fatherPanel) id <- fatherNames(ped)[j]
			if(motherPanel) id <- motherNames(ped)[j]
			if(offsprPanel) id <- sampleNames(ped)[j]
			lrr.segments <- lrr.segments[sampleNames(lrr.segments) %in% id & as.character(chromosome(lrr.segments)) == as.character(chromosome(range)), ]
			if(length(lrr.segments)>0){
				lsegments(x0=start(lrr.segments)/1e6,
					  x1=end(lrr.segments)/1e6,
					  y0=elementMetadata(lrr.segments)$seg.mean,
					  y1=elementMetadata(lrr.segments)$seg.mean, lwd=2, col=segment.col)
			}
		}
	}
	if(panel.number() > 1){ ## label axis for BAFs
		at <- c(blim[1], mean(c(blim[2], blim[1])), blim[2])
		panel.axis("right", at=at, labels=c(0, 0.5, 1), text.col="blue", line.col="blue", half=FALSE, text.cex=0.7)
	}
}



