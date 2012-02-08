xyplotTrioLrrBaf <- function(rd, object, frame, lrr.segments, md.segments, ...){
	index <- seq_len(nrow(rd))
	df <- foreach(i=index, .combine="rbind") %do% {
		dataFrameFromRange2(range=rd[i, ],
				    object=object,
				    frame=frame,
				    range.index=i)
	}
	df$range <- factor(paste("range", df$range), ordered=TRUE, levels=unique(paste("range", df$range)))
	index <- split(seq_len(nrow(df)), df$range)
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
	return(figs)
}


xyplotTrioListLrrBaf <- function(rd, md, object, frame,
				 lrr.segments, md.segments, ...){
	## assume rd is one range
	object <- object[[chromosome(rd)]]
	marker.index <- subjectHits(findOverlaps(rd, featureData(object), maxgap=frame))
	trio.index <- match(sampleNames(rd), sampleNames(object))
	object <- object[marker.index, trio.index]
	md <- md[[chromosome(rd)]]
	md <- md[marker.index, trio.index, drop=FALSE]
	mindist(object) <- md
 	xyplotTrioLrrBaf(rd=rd, object=object, frame=frame, lrr.segments=lrr.segments, md.segments=md.segments, ...)
}

xypanelTrioBaf <- function(x, y,
			   memberId,
			   baf,
			   is.snp,
			   range,
			   lrr.segments,
			   md.segments,
			   col.hom="grey20",
			   fill.hom="lightblue",
			   col.het="grey20" ,
			   fill.het="salmon",
			   col.np="grey20",
			   fill.np="grey60",
			   show.state=TRUE,
			   cex.state=1,
			   col.state="blue",
			   cex.pch=0.3,
			   ped,
			   ..., subscripts){
	panel.grid(v=0, h=4, "grey", lty=2)
	panel.xyplot(x[1], y[1], col="white", cex=cex.pch, ...) ## set it up, but don't plot
	is.snp <- is.snp[subscripts]
	ylim <- current.panel.limits()$ylim
	y[y>ylim[2]] <- ylim[2]

	lpoints(x[!is.snp], y[!is.snp], col=col.np,
		fill=fill.np, cex=cex.pch, ...)
	## use whatever col.hom to color SNPs
	lpoints(x[is.snp], y[is.snp], col=col.hom,
		fill=fill.hom, cex=cex.pch, ...)
	j <- panel.number()
	st <- start(range)[j]/1e6
	lrect(xleft=st, xright=end(range)[j]/1e6,
	      ybottom=-10, ytop=10, ...)
	if(show.state){
		## left justify the label to the start of the range
		y.max <- ylim[2]
		ltext(st, y.max, labels=paste("state", state(range)[j]),
		      adj=c(0,1), cex=cex.state, col=col.state)
	}
	b <- baf[subscripts]
	b[b==2] <- NA
	blim <- c(ylim[1]+0.1, ylim[1]+1.5)
	bnew <- rescale(b, blim[1], blim[2])
	lpoints(x[is.snp], bnew[is.snp], cex=cex.pch, col="blue", ...)
	memberId <- unique(memberId[subscripts])
	##sns <- unique(sampleNames[subscripts])
	if(memberId == "min dist"){
		md.segments <- md.segments[sampleNames(md.segments) %in% sampleNames(range) & chromosome(md.segments) == chromosome(range), ]
		lsegments(x0=start(md.segments)/1e6,
			 x1=end(md.segments)/1e6,
			 y0=mean(md.segments),
			 y1=mean(md.segments), lwd=2, col="grey50")
	} else {
		## range is labeled by offspring id.
		j <- match(sampleNames(range), sampleNames(ped))
		if(memberId == "father"){
			id <- fatherNames(ped)[j]
		} else {
			if(memberId == "mother"){
				id <- motherNames(ped)[j]
			} else {
				id <- sampleNames(range)
			}
		}
		lrr.segments <- lrr.segments[sampleNames(lrr.segments) %in% id & chromosome(lrr.segments) == chromosome(range), ]
		if(nrow(lrr.segments)>0){
			lsegments(x0=start(lrr.segments)/1e6,
				  x1=end(lrr.segments)/1e6,
				  y0=mean(lrr.segments),
				  y1=mean(lrr.segments), lwd=2, col="grey50")
		}
	}
	if(panel.number() > 1){ ## label axis for BAFs
		at <- c(blim[1], mean(c(blim[2], blim[1])), blim[2])
		panel.axis("right", at=at, labels=c(0, 0.5, 1), text.col="blue", line.col="blue", half=FALSE)
	}
}
