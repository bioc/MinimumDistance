catFun2 <- function(rd.query, rd.subject, ...){
	##stopifnot(nrow(rd.query) == nrow(rd.subject)) ## must compare same list size
	ir.q <- IRanges(start(rd.query), end(rd.query))
	ir.s <- IRanges(start(rd.subject), end(rd.subject))
	mm <- findOverlaps(ir.q, ir.s, ...)
	query.index <- queryHits(mm)
	if(length(query.index) == 0){
		return(0)
	}
	subject.index <- subjectHits(mm)
	index <- which(chromosome(rd.query)[query.index] == chromosome(rd.subject)[subject.index] &
		       sampleNames(rd.query)[query.index] == sampleNames(rd.subject)[subject.index])
	if(length(index) > 0){
		query.index <- unique(query.index[index])
		p <- length(query.index)/nrow(rd.query)
		if(p > 1) browser()
	} else p <- 0
	return(p)
}

splitByDistance <- function(x, thr=100e3){
	d <- diff(x)
	if(all(d < thr)) return(rep(0, length(x)))
	f <- c(0, cumsum(d > thr))
	tab.f <- table(f)
	## combine regions if number of markers is very small
	while(any(tab.f < 1000) & length(x) > 1000){
		j <- which(tab.f < 1000)[[1]]
		factor.val <- as.integer(names(tab.f)[j])
		if(factor.val < max(f)){
			f[f==factor.val] <- factor.val+1
		} else {
			f[f==factor.val] <- factor.val-1
		}
		tab.f <- table(f)
	}
	return(f)
}

splitIndicesByLength2 <- function(x, MIN.LENGTH=1000, ...){
	f <- splitIndicesByLength(x, ...)
	l <- sapply(f, length)
	L <- length(l)
	l <- l[L]
	if(l < MIN.LENGTH & length(f) > 1){
		f[[L-1]] <- c(f[[L-1]], f[[L]])
		f <- f[-L]
	}
	return(f)
}



discAtTop <- function(ranges.query, ranges.subject, verbose=TRUE,...){
	ir.q <- IRanges(start(ranges.query), end(ranges.query))
	ir.s <- IRanges(start(ranges.subject), end(ranges.subject))
	mm <- findOverlaps(ir.q, ir.s,...)
	query.index <- queryHits(mm)
	subject.index <- subjectHits(mm)
	index <- which(chromosome(ranges.query)[query.index] == chromosome(ranges.subject)[subject.index] &
		       sampleNames(ranges.query)[query.index] == sampleNames(ranges.subject)[subject.index])
	query.index <- unique(query.index[index])
	##subject.index <- unique(subject.index[index])
	notOverlapping.index <- seq(length=nrow(ranges.query))[!seq(length=nrow(ranges.query)) %in% query.index]
	res <- ranges.query[notOverlapping.index, ]
	return(res)
}

concAtTop <- function(ranges.query, ranges.subject, list.size, verbose=TRUE, ...){
	p <- rep(NA, length(list.size))
	pAny1 <- rep(NA, length(list.size))
	pAny2 <- rep(NA, length(list.size))
	if(verbose) {
		message("Calculating the proportion of ranges in common for the first ", max(list.size), " ranges")
		pb <- txtProgressBar(min=0, max=length(p), style=3)
	}
	for(i in seq_along(list.size)){
		if(verbose) setTxtProgressBar(pb, i)
		L <- list.size[i]
		p[i] <- catFun2(ranges.query[seq(length=L), ], ranges.subject[seq(length=L), ], ...)
		pAny1[i] <- catFun2(ranges.query[seq(length=L), ], ranges.subject, ...)
		pAny2[i] <- catFun2(ranges.subject[seq(length=L), ], ranges.query, ...)
	}
	if(verbose) close(pb)
	res <- list(p=p, pAny.queryList=pAny1, pAny.subjectList=pAny2)
	names(res) <- c("cat", "topMD", "topPenn")
	return(res)
}


correspondingCall <- function(ranges.query, ranges.subject, subject.method){
	overlap <- findOverlaps(ranges.query, ranges.subject)
	subj.index <- subjectHits(overlap)
	quer.index <- queryHits(overlap)
	## what are the chromosomes for the subject hits
	index <- which(chromosome(ranges.query)[quer.index] == chromosome(ranges.subject)[subj.index] &
		       sampleNames(ranges.query)[quer.index] == sampleNames(ranges.subject)[subj.index])
	if(length(index) == 0) return("no overlap")
	matching.index <- subj.index[index]
	res <- ranges.subject[matching.index, ]
	if(!missing(subject.method)) res$method <- subject.method
	return(res)
}

isDeletion <- function(x){
	if(length(grep("-", x)) > 0){
		tmp <- strsplit(x, "_")[[1]]
		state <- substr(tmp, 3, 3)
		state <- ifelse(any(state < 3), TRUE, FALSE)
	} else{
		state <- as.integer(substr(x, 3, 3))
		state <- ifelse(state < 3, TRUE, FALSE)
	}
	state
}

overlapsCentromere <- function(myranges){
	##require(SNPchip)
	data(chromosomeAnnotation, package="SNPchip")
	chromosomeAnnotation <- get("chromosomeAnnotation")
	centromere.ranges <- RangedData(IRanges(chromosomeAnnotation[, "centromereStart"],
						chromosomeAnnotation[, "centromereEnd"]),
					chrom=rownames(chromosomeAnnotation))
	myranges.bak <- myranges
	chrom <- unique(myranges$chrom)
	overlaps.centromere <- rep(NA, nrow(myranges))
	for(CHR in chrom){
		centromere.ir <- IRanges(start(centromere.ranges)[CHR],
					 end(centromere.ranges)[CHR])
		ix <- which(myranges$chrom==CHR)
		ir <- IRanges(start(myranges)[ix],
			      end(myranges)[ix])
		overlaps.centromere[ix] <- countOverlaps(ir, centromere.ir) > 0
	}
	return(overlaps.centromere)
}

getRefGene <- function(filename="~/Data/Downloads/hg18_refGene.txt"){
	colClasses <- c("integer", "character", "character", "factor",
			"integer", "integer",
			"integer", "integer",
			"integer",
			"character", "character",
			"integer", rep("character", 4))
	tmp <- read.delim(filename, header=FALSE,
			  colClasses=colClasses)
	tmp <- tmp[, c(2:6, 13)]
	colnames(tmp) <- c("NM", "chrom", "strand", "start", "end", "gene_name")
	chrom <- sapply(tmp$chrom, function(x) strsplit(x, "chr")[[1]][2])
	tmp$chrom <- chromosome2integer(chrom)
	tmp <- tmp[!is.na(tmp$chrom), ]
	refGene <- RangedData(IRanges(tmp$start, tmp$end),
			      chrom=tmp$chrom,
			      strand=tmp$strand,
			      NM=tmp$NM,
			      gene_name=tmp$gene_name)
	refGene
}

combineRanges <- function(deletion.ranges, amp.ranges){
	state <- deletion.ranges$state
	hemizygous.states <- c("332", "432", "342")
	homozygous.states <- c("331", "321", "231", "431", "341", "441", "221")
	deletion.ranges <- deletion.ranges[state %in% hemizygous.states | state %in% homozygous.states, ]
	amp.ranges <- amp.ranges[, colnames(amp.ranges) %in% colnames(deletion.ranges)]
	index <- match(colnames(amp.ranges), colnames(deletion.ranges))
	deletion.ranges2 <- deletion.ranges[,  index]
	stopifnot(all.equal(colnames(deletion.ranges2), colnames(amp.ranges)))
	ranges.all <- RangedData(IRanges(c(start(deletion.ranges2), start(amp.ranges)),
					 c(end(deletion.ranges2), end(amp.ranges))),
				 id=c(deletion.ranges2$id, amp.ranges$id),
				 chrom=c(deletion.ranges2$chrom, amp.ranges$chrom),
				 num.mark=c(deletion.ranges2$num.mark, amp.ranges$num.mark),
				 seg.mean=c(deletion.ranges2$seg.mean, amp.ranges$seg.mean),
				 state=c(deletion.ranges2$state, amp.ranges$state))
	ranges.all
}


pruneByFactor <- function(range.object, f, verbose=FALSE){
	rd <- list()
	id.chr <- paste(sampleNames(range.object), chromosome(range.object), sep="_")
	ff <- unique(id.chr)
	##for(i in seq_along(unique(range.object$id))){
	if(verbose){
		message("Pruning ", length(ff), " files.")
		pb <- txtProgressBar(min=0, max=length(ff), style=3)
	}
	## do for each element in a GRangesList
	for(i in seq_along(ff)){
		if(verbose) setTxtProgressBar(pb, i)
		##id <- unique(range.object$id)[i]
		##(index <- which(range.object$id == id))
		index <- which(id.chr==ff[i])
		##trace(combineRangesByFactor, browser)
		rd[[i]] <- combineRangesByFactor(range.object[index, ], f=f[index])
	}
	if(verbose) close(pb)
	##ok <- tryCatch(stack(RangedDataList(rd)), error=function(e) FALSE)
	##rd <- GRangesList(rd)
	rd <- unlist(GRangesList(rd))
	##rd <- stackRangedDataList(rd)
##	names(rd) <- unique(sampleNames(range.object))
##	if(!is(ok, "RangedData")) {
##		message("trouble combining RangedData objects.  Returning list")
##		ok <- rd
##	} else {
##		j <- match("sample",colnames(ok))
##		if(length(j) == 1)
##			ok <- ok[, -j]
##	}
	return(rd)
}

combineRangesByFactor <- function(range.object, f){
	##range.object <- range.object[!is.na(state(range.object)), ]
	i <- which(is.na(f))
	j <- 1
	while(length(i) > 0){
		if(is.na(f[1])){
			f[1] <- f[2]
		} else {
			f[is.na(f)] <- f[i-1]
		}
		i <- which(is.na(f))
		j <- j+1
		if(j > 10) stop("too many na's in f")
	}
	##stopifnot(all(!is.na(f)))
	ff <- cumsum(c(0, abs(diff(as.integer(as.factor(f))))))
	if(!any(duplicated(ff))) {
		return(range.object)
	}
	for(i in seq_along(unique(ff))){
		x <- unique(ff)[i]
		if(sum(ff==x) == 1) next()
		index <- which(ff==x)
		min.index <- min(index)
		max.index <- max(index)
		end(range.object)[index] <- max(end(range.object)[index])
		emd <- elementMetadata(range.object)
		emd$lik.state[index] <- sum(emd$lik.state[index], na.rm=TRUE)
		emd$seg.mean[index] <- sum((numberProbes(range.object)[index]*emd$seg.mean[index]), na.rm=TRUE)/sum(numberProbes(range.object)[index], na.rm=TRUE)
		emd$numberProbes[index] <- sum(numberProbes(range.object)[index], na.rm=TRUE)
		emd$lik.norm[index] <- sum(emd$lik.norm[index], na.rm=TRUE)
		elementMetadata(range.object) <- emd
		j <- seq_len(length(range.object))
		index <- index[-1]
		j <- j[-index]
		if(length(j) == 0){
			stop()
		}
		ff <- ff[j]
		range.object <- range.object[j, ]
	}
	return(range.object)
}

madVsCoverage <- function(lambda=0.1, MIN=1, MAX=4, coverage=3:100){
	p <- lambda*exp(-lambda*coverage) ## 0 - 0.04 (Pr (X=x)
	b <- 1/(MAX - MIN)
	a <- MIN * b
	numberMads <- ((p-min(p))/(max(p)-min(p)) + a)/b
	list(x=coverage, y=numberMads)
}

thresholdSegMeans <- function(ranges.object, ylim){
	ranges.object$seg.mean[ranges.object$seg.mean < ylim[1]] <- ylim[1]
	ranges.object$seg.mean[ranges.object$seg.mean > ylim[2]] <- ylim[2]
	ranges.object
}

combine.data.frames <- function(dist.df, penn.df){
	if(is.null(dist.df) & is.null(penn.df)) return(NULL)
	if(is.null(dist.df)) dist.df <- penn.df[integer(0), ]
	if(is.null(penn.df)) penn.df <- dist.df[integer(0), ]
	combined.df <- rbind(dist.df, penn.df)
	combined.df <- combined.df[order(combined.df$chr), ]
	return(combined.df)
}



##offspring.hemizygousPenn <- function() c("332", "432", "342", "442")
offspring.hemizygousPenn <- function(){
	tmp <- expand.grid(c(1,3,5,6), c(1,3,5,6), 1)
	dels <- paste(tmp$Var1, tmp$Var2, tmp$Var3, sep="")
	dels <- dels[-1]
}
##offspring.hemizygous <- function() c("221", "321", "231", "441", "341", "431")
offspring.hemizygous <- function() {
	tmp <- expand.grid(c(0,2,3,4), c(0,2,3,4), 1)
	dels <- paste(tmp$Var1, tmp$Var2, tmp$Var3, sep="")
	dels <- dels[-1]
	dels
}
offspring.homozygous <- function(){
	tmp <- expand.grid(c(1,2,3,4), c(1,2,3,4), 0)
	dels <- paste(tmp$Var1, tmp$Var2, tmp$Var3, sep="")
	dels
}
##offspring.homozygous <- function() c("220", "210", "120", "320", "230", "330", "110", "310")
deletionStates <- function(){
	st1 <- offspring.hemizygous()
	st2 <- offspring.homozygous()
	as.integer(c(st1,st2))
}
##duplicationStates <- function() as.integer(c("224", "223", "113", "114", "013", "014", "103", "104", "213", "214", "123", "124", "013", "014", "103", "104", "203", "204", "023", "024"))
duplicationStates <- function(){
	tmp <- expand.grid(c(0,1,2,4), c(0,1,2,4), 3)
	sdups <- paste(tmp$Var1, tmp$Var2, tmp$Var3, sep="")
	sdups <- sdups[-1]
	tmp <- expand.grid(c(0,1,2,3), c(0,1,2,3), 4)
	ddups <- paste(tmp$Var1, tmp$Var2, tmp$Var3, sep="")
	ddups <- ddups[-1]
	c(sdups, ddups)
}
##duplicationStatesPenn <- function() as.integer(c("335", "336", "225", "226", "115", "116", "125", "215", "325", "235", "125", "215", "315", "135"))
duplicationStatesPenn <- function() {
	tmp <- expand.grid(c(1,2,3,6), c(1,2,3,6), 5)
	sdups <- paste(tmp$Var1, tmp$Var2, tmp$Var3, sep="")
	sdups <- sdups[-1]
	tmp <- expand.grid(c(1,2,3,6), c(1,2,3,6), 5)
	ddups <- paste(tmp$Var1, tmp$Var2, tmp$Var3, sep="")
	ddups <- ddups[-1]
	c(sdups, ddups)
##	dups.penn <- expand.grid(c(1,2,3,5,6), c(1,2,3,5,6), c(5,6))
##	paste(dups.penn$Var1, dups.penn$Var2, dups.penn$Var3)
}
isDenovo <- function(states) states %in% c(duplicationStates(), deletionStates())

##offspring.hemizygous <- function() c("332", "532", "352", "552")
##offspring.homozygous <- function() c("331", "321", "231", "531", "351", "551", "221", "521")
##duplicationStates <- function() as.integer(c("336", "335", "225", "226", "116", "115", "125", "126", "215", "216", "325", "326", "235", "236", "125", "126", "215", "216", "315", "316", "135", "136"))
####duplicationStatesPenn <- function() as.integer(c("335", "225", "115", "125", "215", "325", "235", "125", "215", "315", "135"))
##isDenovo <- function(states) states %in% c(duplicationStates(), deletionStates())

calculateChangeSd <- function(coverage=1:500, lambda=0.05, a=0.2, b=0.025)
	a + lambda*exp(-lambda*coverage)/b

pruneMD <- function(genomdat,
		  range.object,
		  physical.pos,
		  ##trimmed.SD, ##
		  lambda=0.05,
		  MIN.CHANGE=0.1,
		  SCALE.EXP=0.02,
		  MIN.COVERAGE=3,
		  weighted=FALSE,
		  weights=NULL) {
	if(length(unique(range.object$id)) != 1) stop("multiple ids in range.object")
	if(length(unique(chromosome(range.object))) > 1) stop("Multiple chromosomes in range.object")
	##change.SD <- trimmed.SD*change.SD
	genomdat <- as.numeric(genomdat)/100
	coverage <- range.object$num.mark
	trimmed.SD <- max(mad(genomdat, na.rm=TRUE), .15)
	##trimmed.SD <- unique(range.object$mindist.mad)
	##stopifnot(length(trimmed.SD)==1)
	coverage <- coverage[-length(coverage)]
	if(FALSE){
		numberSds <- calculateChangeSd(coverage=3:100, lambda=lambda, a=MIN.CHANGE, b=SCALE.EXP)
		y <- MIN.CHANGE+lambda*exp(-lambda*(3:100))/SCALE.EXP
		plot(3:100, y, ylab="number of MADs", xlab="coverage")
	}
	##thrSD <- calculateChangeSd(coverage, lambda, trimmed.SD, change.SD)
		##change.SD <- change.SD  ##Thresholds for right cutpoint
	##cpt.loc <- cumsum(lseg) ## indices of the cutpoints same as coverage.
	cpt.loc <- range.object$end.index
	sdundo <- TRUE
	while(sdundo) {
		k <- length(cpt.loc)
		if (k>1) {
			coverage <- diff(c(0, cpt.loc))
			coverage <- coverage[-length(coverage)]
			##
			##  number of sds as a function of coverage
			##  -- segments with high coverage have small y
			##
			requiredNumberSd <- calculateChangeSd(coverage=coverage, lambda=lambda, a=MIN.CHANGE, b=SCALE.EXP)
			##
			## number of standard deviations
			segments0 <- cbind(c(1,1+cpt.loc[-k]),cpt.loc)
			## median copy number for each segment
			segmed <- apply(segments0, 1, function(i,x) {median(x[i[1]:i[2]], na.rm=T)}, genomdat)
			## absolute copy number difference of adjacent segments
 			##adsegmed <- abs(diff(segmed))
			adsegmed <- abs(diff(segmed))
			## number of standard deviations of observed shift
			empiricalNumberSd <- adsegmed/trimmed.SD
			if(any(empiricalNumberSd < requiredNumberSd | coverage < MIN.COVERAGE)){
				## drop order: coverage then distance
				##i <- which(adsegmed < thrSD | coverage < MIN.COVERAGE)
				i <- which(empiricalNumberSd < requiredNumberSd | coverage < MIN.COVERAGE)
				if(length(i) > 1){
					i <- i[order(coverage[i], adsegmed[i], decreasing=FALSE)[1]]
				}
				cpt.loc <- cpt.loc[-i]
			} else {
				sdundo <- FALSE
			}
		} else {
			sdundo <- FALSE
		}
	}
	lseg <- diff(c(0,cpt.loc)) ## back to coverage
	## update segment means
	segmeans <- 0*lseg
	ll <- uu <- 0
	for(i in 1:length(lseg)) {
		uu <- uu + lseg[i]
		if (weighted) {
			segmeans[i] <- sum(genomdat[(ll+1):uu]*weights[(ll+1):uu])/sum(weights[(ll+1):uu])
		} else {
			segmeans[i] <- mean(genomdat[(ll+1):uu], na.rm=TRUE)
		}
		ll <- uu
	}
	segments0 <- cbind(c(1,1+cpt.loc[-k]),cpt.loc)
	starts <- physical.pos[segments0[, 1]]
	ends <- physical.pos[segments0[, 2]]
	if(length(ends) < length(starts)) ends <- c(ends, max(end(range.object)))
	id <- unique(range.object$id)
	res <- RangedDataCBS(IRanges(starts, ends),
			     sampleId=id,
			     chrom=unique(range.object$chrom),
			     coverage=lseg,
			     seg.mean=segmeans,
			     mindist.mad=mad(genomdat, na.rm=TRUE))
	return(res)
}

## pdf of standard normal
## the msm package has this stuff, but it seemed slow...
phi <- function(x, mu, sigma) dnorm(x, mu, sigma)
## cdf of standard normal
Phi <- function(x, mu, sigma) pnorm(x, mu, sigma)
## pdf of truncated normal on support [0, 1]
tnorm <- function(x, mean, sd, lower=0, upper=1){
	res <- phi(x, mean, sd)/(Phi(upper, mean, sd)-Phi(lower, mean, sd))
	ind <- which(x < lower | x > upper)
	if(any(ind)){
		res[ind] <- 0
	}
	res
}
TN <- tnorm


addRangeIndex <- function(id, trioSet, ranges){
	ranges <- ranges[sampleNames(ranges) %in% id, ]
	stopifnot(nrow(ranges) > 0)
	stopifnot(id %in% sampleNames(trioSet))
	ir1 <- IRanges(start=position(trioSet), end=position(trioSet))
	ir2 <- IRanges(start(ranges), end(ranges))
	mm <- findOverlaps(ir1, ir2)
	## there should be no query that is in more than 1 subject
	qhits <- queryHits(mm)
	shits <- subjectHits(mm)
	right.chromosome <- chromosome(ranges)[shits] == chromosome(trioSet)[qhits]
	qhits <- qhits[right.chromosome]
	shits <- shits[right.chromosome]
	range.index <- rep(NA, nrow(trioSet))
	##fData(object)$range.index <- NA
	##fData(object)$range.index[qhits] <- shits
	range.index[qhits] <- shits
	if(sum(table(range.index)) != nrow(trioSet)){
		message("# of markers in the ranges not equal to total number of markers")
		browser()
	}
	return(range.index)
}

pHet <- function(i, id, trioSet){
	j <- match(id, sampleNames(trioSet))
	stopifnot(length(j) > 0)
	is.ff <- is(baf(trioSet), "ff")
	if(is.ff){
		open(baf(trioSet))
	}
	b <- baf(trioSet)[i, j, 3]
	if(is.ff){
		close(baf(trioSet))
	}
	mean(b > 0.4 & b < 0.6, na.rm=TRUE)
}
meanLogR <- function(i, id, trioSet){
	j <- match(id, sampleNames(trioSet))
	stopifnot(length(j) > 0)
	is.ff <- is(lrr(trioSet), "ff")
	if(is.ff){
		open(lrr(trioSet))
	}
	r <- lrr(trioSet)[i, j, 3]
	if(is.ff){
		close(lrr(trioSet))
	}
	mean(r, na.rm=TRUE)
}


LikSet <- function(trioSet, pedigreeData, id, CHR, ranges){
	is.ff <- is(lrr(trioSet), "ff")
	if(missing(id)) id <- sampleNames(trioSet)[1]
	if(is.ff){
		open(baf(trioSet))
		open(lrr(trioSet))
	}
	i <- match(id, sampleNames(trioSet))
	stopifnot(length(i) == 1)
	## the trios are in the same order as the sampleNames of the trioSetList object
	## validity methods for the class ensure that this is correct
	indNames <- as.character(trios(pedigreeData)[i,])
	##offspring.id <- id[id %in% offspringNames(trioSet)]
	##i <- match(offspring.id, sampleNames(trioSet))
	##i <- match(id[["O"]], offspringNames(trioSet))
	mads <- mad(trioSet)[i, ]
	##S <- length(states)
	loglik <- array(NA, dim=c(2, nrow(trioSet), 3, 5))
	dimnames(loglik) <- list(c("logR", "baf"),
				 featureNames(trioSet),
				 indNames,
				 0:4)
	object <- new("LikSet",
		      logR=as.matrix(lrr(trioSet)[ ,i,]),
		      BAF=as.matrix(baf(trioSet)[ ,i , ]),
		      featureData=featureData(trioSet),
		      loglik=loglik)
	object$MAD <- mads
	fData(object)$range.index <- NA
	##tmp=findOverlaps(featureData(object), ranges)
	fo <- findOverlaps(ranges, featureData(object))
	i1 <- subjectHits(fo)
	i2 <- queryHits(fo)
	fData(object)$range.index[i1] <- i2
	if(any(is.na(range.index(object)))){
		msg <- paste("Segmentation was run on chunks of the data for which the markers are less than 75kb apart.\n",
			     "When log R ratios are missing at the boundaries of the partioned data, not all markers \n",
			     "will be covered by a segment.\n")
		if(is.null(.GlobalEnv[[".warningMessageAlreadyDisplayed"]])){
			warning(msg)
			.GlobalEnv[[".warningMessageAlreadyDisplayed"]] <- TRUE
		}
		object <- object[!is.na(range.index(object)), ]
	}
	## NA values occur if there are ranges that do not overlap the
	## marker locations in object.
##	ir1 <- IRanges(start=position(object), end=position(object))
##	ir2 <- IRanges(start(ranges), end(ranges))
##	mm <- findOverlaps(ir1, ir2)
##	subject.index <- subjectHits(mm)
##	## there should be no query that is in more than 1 subject
##	qhits <- queryHits(mm)
##	shits <- subjectHits(mm)
##	fData(object)$range.index <- NA
##	fData(object)$range.index[qhits] <- shits
	if(is.ff){
		close(baf(trioSet))
		close(lrr(trioSet))
	}
	return(object)
}

fillInMissing <- function(rangeIndex){
	if(!any(is.na(rangeIndex))) return(rangeIndex)
	if(sum(is.na(rangeIndex)) > 1000 & length(unique(rangeIndex[!is.na(rangeIndex)])) == 1){
		## for calculating the posterior for a single range
		## essentially, ignoring what comes before and what comes after
		ii <- range(which(!is.na(rangeIndex)))
		if(ii[1] > 1){
			rangeIndex[1:(ii[1]-1)] <- 0
		}
		if(ii[2] < length(rangeIndex)){
			rangeIndex[(ii[2]+1):length(rangeIndex)] <- 2
		}
	} else {
		ii <- which(is.na(rangeIndex))
		if(max(ii) < length(rangeIndex)){
			rangeIndex[ii] <- rangeIndex[ii+1]
		} else{
			rangeIndex[max(ii)] <- rangeIndex[max(ii)-1]
			if(any(ii != max(ii))){
				iii <- ii[ii != max(ii)]
				rangeIndex[iii] <- rangeIndex[iii+1]
			}
		}
	}
	return(rangeIndex)
}

rowMAD <- function(x, y, ...){
	##notna <- !is.na(x)
	mad <- 1.4826*rowMedians(abs(x-rowMedians(x, ...)), ...)
	return(mad)
}

dups.penn <- expand.grid(c(1,2,3,5,6), c(1,2,3,5,6), c(5,6))

trioStates <- function(states=0:4){
	trio.states <- as.matrix(expand.grid(states, states, states))
	index <- which(trio.states[, 1] == 0 & trio.states[, 2] == 0 & trio.states[, 3] > 0)
	##trio.states <- trio.states+1
	colnames(trio.states) <- c("F", "M", "O")
	## 125 possible
	## remove 00 > 0 as possibilities
	trio.states <- trio.states[-index, ]
}

trioStateNames <- function(trio.states){
	paste(paste(trio.states[,1], trio.states[,2], sep=""), trio.states[,3], sep="")
}


transitionProbability <- function(nstates=5, epsilon=1-0.999){
	off.diag <- epsilon/(nstates-1)
	tpm <- matrix(off.diag, nstates, nstates)
	diag(tpm) <- 1-epsilon
	tpm
}



readTable1 <- function(states=0:4, a=0.0009){
	S <- length(states)
	tmp <- array(NA, dim=rep(S,3))
	dimnames(tmp) <- list(paste("F", states, sep=""),
			      paste("M", states, sep=""),
			      paste("O", states, sep=""))
	tmp["F0", "M0", ] <- c(1, rep(0,4))
	tmp["F0", "M1", ] <- c(0.5, 0.5, 0, 0, 0)
	tmp["F0", "M2", ] <- c(0.5*a, 1-a, 0.5*a, 0, 0)
	tmp["F0", "M3", ] <- c(0.5*a, 0.5*(1-a), 0.5*(1-a), 0.5*a, 0)
	tmp["F0", "M4", ] <- c(0, 0.25, 0.5, 0.25, 0)

	tmp["F1", "M1", ] <- c(0.25, 0.5, 0.25, 0, 0)
	tmp["F1", "M2", ] <- c(0.25, 0.5-0.25*a, 0.5-0.25*a, 0.25*a, 0)
	tmp["F1", "M3", ] <- c(0.25*a, 0.25, 0.5*(1-a), 0.25, 0.25*a)
	tmp["F1", "M4", ] <- c(0, 0.125, 0.375, 0.375, 0.125)

	tmp["F2", "M2", ] <- c(0.25*a^2, a*(1-a), (1-a)^2 + 0.5*a^2, a*(1-a), 0.25*a^2)
	tmp["F2", "M3", ] <- c(0.25*a^2, 0.75*a*(1-a), 0.5*(1-a)^2+0.25*a*(1-a)+0.25*a^2,
			       0.5*(1-a)^2+0.25*a*(1-a)+0.25*a^2,
			       0.75*a*(1-a)+0.25*a^2)
	tmp["F2", "M4", ] <- c(0,0.125*a, 0.25, 0.5-0.25*a,0.25+0.125*a)

	tmp["F3", "M3", ] <- c(0.25^2, 0.5*a*(1-a), 0.5*a*(1-a)+0.25*(1-a)^2, 0.5*(1-a)^2+0.5^2,
			       0.25*(1-a)^2 + a*(1-a)+0.25^2)
	tmp["F3", "M4", ] <- c(0, 0.125*a, 0.125*(1+a), 0.125*(3-2*a), 0.5)

	tmp["F4", "M4", ] <- c(0, 0, 0.0625, 0.25, 0.6875)
	return(tmp)
}

lookUpTable1 <- function(state, table1){
  index <- state+1
  p <- table1[index[1], index[2], index[3]]
  if(is.na(p)){
    p <- table1[index[2], index[1], index[3]]
  }
  p
}

vectorizeTable1 <- function(table1, stateMatrix){
  setNames(apply(stateMatrix, 1, lookUpTable1, table1=table1), rownames(stateMatrix))
}

lookUpTable3 <- function(table3, state.prev, state.curr){
	f1 <- state.prev[1]
	f2 <- state.curr[1]
	m1 <- state.prev[2]
	m2 <- state.curr[2]
	o1 <- state.prev[3]
	o2 <- state.curr[3]
	return(table3[f1, f2, m1, m2, o1, o2])
}

jointProb <- function(state,
                      param,
		      state.prev,
		      log.lik){
  if(is.null(state.prev)){
    result <- setNames(loglik222(param, LLT=log.lik), "222")
    return(result)
  }
  ##browser()
  state_name <- state
  state_index <- state(param)[state_name, ] + 1L
  prev_index <- state.prev + 1L
  log.pi <- log(initialStateProb(param))
  tau <- transitionProb(param)
  pi.f <- exp(log.pi[state_index[1]])
  pi.m <- exp(log.pi[state_index[2]])
  tau.o <- tau[prev_index[3], state_index[3]]
  tau.m <- tau[prev_index[2], state_index[2]]
  tau.f <- tau[prev_index[1], state_index[1]]
  p.00 <- lookUpTable3(table3(param), prev_index, state.curr=state_index) ## both Mendelian
  p.10 <- 1/5*table1(param)[state_name]
  ##p.10 <- 1/5*lookUpTable1(table1, state) ## previous non-Mendelian * current Mendelian
  prev_name <- paste0(state.prev, collapse="")
  p.01 <- table1(param)[prev_name]*1/5 ## previous Mendlian * current non-mendelian
  ##p.01 <- lookUpTable1(table1, state.prev) * 1/5 ## previous Mendlian * current non-mendelian
  p.11 <- 1/5*tau.o
  prNM <- prNonMendelian(param)
  piI0 <- 1-prNM ##prob.nonMendelian
  p.NM <- piI1 <- prNM
  ## setting this to a small value will favor '2,2,1' versus '3,3,1' (for example)
  ## setting prob.nonMendelian smaller would not have an effect
  tauI.11 <- tauI.00 <- 1-.01
  tauI.10 <- tauI.01 <- 1-tauI.11
  pr.off <- p.NM*(p.11*tauI.11 + p.10*tauI.10) + (1-p.NM)*(p.00*tauI.00+p.01*tauI.01)
  LL <- diag(log.lik[, state_index])
  log.lik <- setNames(sum(LL) + log(pr.off*tau.m*pi.m*tau.f*pi.f),
                      state_name)
  log.lik
}

xypanelMD <- function(x, y,
		      id,
		      gt,
		      is.snp,
		      range,
		      cex,
		      col.hom="grey20",
		      fill.hom="lightblue",
		      col.het="grey20" ,
		      fill.het="salmon",
		      col.np="grey20",
		      fill.np="grey60",
		      show.state=TRUE,
		      lrr.segs,
		      md.segs,
		      ..., subscripts){
	xypanel(x, y,
                gt,
                is.snp,
                range,
                col.hom=col.hom,
                fill.hom=fill.hom,
                col.het=col.het,
                fill.het=fill.het,
                col.np=col.np,
                fill.np=fill.np,
                show.state, cex=cex, ..., subscripts=subscripts)
	id <- unique(id[subscripts])
	range <- range[1, ]
	CHR <- chromosome(range)
	##stopifnot(length(CHR)==1)
	if(id != "min dist" & !missing(lrr.segs)){
		cbs.sub <- lrr.segs[sampleNames(lrr.segs)==as.character(id) & chromosome(lrr.segs)==CHR, ]
		segments <- TRUE && nrow(cbs.sub) > 0
	} else segments <- FALSE
	if(!missing(md.segs) & id == "min dist"){
		cbs.sub <- md.segs[sampleNames(md.segs) %in% sampleNames(range), ]
		cbs.sub <- cbs.sub[chromosome(cbs.sub) == chromosome(range), ]
		##cbs.sub$seg.mean <- -1*cbs.sub$seg.mean
		segments.md <- TRUE && nrow(cbs.sub) > 0
	} else segments.md <- FALSE
	if(segments | segments.md){
		##if(missing(ylimit)) ylimit <- range(y, na.rm=TRUE) ##else ylim <- ylimit
		ylimit <- current.panel.limits()$ylim
		if(nrow(cbs.sub) > 0){
			index <- which(cbs.sub$seg.mean < ylimit[1])
			if(length(index) > 0)
				cbs.sub$seg.mean[index] <- ylimit[1] + 0.2
			index <- which(cbs.sub$seg.mean > ylimit[2])
			if(length(index) > 0)
				cbs.sub$seg.mean[index] <- ylimit[2] - 0.2
			panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2, col="black")#gp=gpar("lwd"=2))
		}
	}
}

xypanelMD2 <- function(x, y,
		       id,
		       gt,
		       is.snp,
		       range,
		       show.state=TRUE,
		       lrr.segs,
		       md.segs,
		       col,
		       cex=1,
		       cex.state=1,
		       col.state="blue",
		       ..., subscripts){
	panel.grid(v=0, h=4, "grey", lty=2)
	panel.xyplot(x, y, cex=cex, col=col[subscripts], ...)
	##lpoints(x, y, col=col[subscripts], cex=cex, ...)
##	lpoints(x[!is.snp], y[!is.snp], col=col.np, cex=cex, ...)
	id <- unique(id[subscripts])
	range <- range[1, ]
	CHR <- chromosome(range)
	##stopifnot(length(CHR)==1)
	if(id != "min dist" & !missing(lrr.segs)){
		cbs.sub <- lrr.segs[sampleNames(lrr.segs)==as.character(id) & chromosome(lrr.segs)==CHR, ]
		segments <- TRUE && nrow(cbs.sub) > 0
	} else segments <- FALSE
	if(!missing(md.segs) & id == "min dist"){
		cbs.sub <- md.segs[sampleNames(md.segs) %in% sampleNames(range), ]
		cbs.sub <- cbs.sub[chromosome(cbs.sub) == chromosome(range), ]
		##cbs.sub$seg.mean <- -1*cbs.sub$seg.mean
		segments.md <- TRUE && nrow(cbs.sub) > 0
	} else segments.md <- FALSE
	if(segments | segments.md){
		##if(missing(ylimit)) ylimit <- range(y, na.rm=TRUE) ##else ylim <- ylimit
		ylimit <- current.panel.limits()$ylim
		if(nrow(cbs.sub) > 0){
			index <- which(cbs.sub$seg.mean < ylimit[1])
			if(length(index) > 0)
				cbs.sub$seg.mean[index] <- ylimit[1] + 0.2
			index <- which(cbs.sub$seg.mean > ylimit[2])
			if(length(index) > 0)
				cbs.sub$seg.mean[index] <- ylimit[2] - 0.2
			panel.segments(x0=start(cbs.sub)/1e6, x1=end(cbs.sub)/1e6, y0=cbs.sub$seg.mean, y1=cbs.sub$seg.mean, lwd=2, col="black")#gp=gpar("lwd"=2))
		}
	}
	j <- panel.number()
	st <- start(range)[j]/1e6
	lrect(xleft=st, xright=end(range)[j]/1e6,
	      ybottom=-10, ytop=10, ...)
	if(show.state){
		## left justify the label to the start of the range
		y.max <- current.panel.limits()$ylim[2]
		ltext(st, y.max, labels=paste("state", state(range)[j]),
		      adj=c(0,1), cex=cex.state, col=col.state)
	}
}

narrow <- function(object, lrr.segs, thr=0.9,
		   mad.minimumdistance, verbose=TRUE,
		   fD, genome) .Defunct("The 'narrow' function is defunct in MinimumDistance. Use narrowRanges instead.")

narrowRanges <- function(object, lrr.segs, thr=0.9,
		   mad.minimumdistance, verbose=TRUE,
		   fD, genome){
  if(missing(fD)) stop("fD not specified. fD must be a list of GenomeAnnotatedDataFrames (if multiple chromosomes are in 'object'), or a single GenomeAnnotatedDataFrame (one chromosome represented in 'object')")
  if(!is(names(mad.minimumdistance), "character")) stop("mad.minimumdistance must be named")
  if(!missing(genome)) metadata(lrr.segs) <- list(genome=genome)
  ix <- match(sampleNames(object), names(mad.minimumdistance))
  object$mindist.mad <- mad.minimumdistance[ix]
  lrr.segs <- lrr.segs[sampleNames(lrr.segs) %in% sampleNames(object), ]
  if(length(unique(chromosome(object))) > 1){
    if(verbose)
      message("narrowing the ranges by chromosome")
    if(!is(fD, "list")) stop("when object contains multiple chromosomes, fD should be a list of GenomeAnnotatedDataFrames")
    chromsInFD <- paste("chr", sapply(fD, function(x) chromosome(x)[1]), sep="")
    indexList <- split(seq_len(length(object)), as.character(chromosome(object)))
    indexList2 <- split(seq_len(length(lrr.segs)), as.character(chromosome(lrr.segs)))
    if(!all.equal(names(indexList), names(indexList2))){
      stop("the chromosomes represented in the minimum distance genomic intervals (object) must be the same as the chromosomes represented in the offspring genomic intervals (lrr.segs)")
    }
    fD <- fD[match(names(indexList), as.character(chromsInFD))]
    if(length(fD) != length(indexList)){
      stop("The list of GenomeAnnotatedDataFrames (argument fD) must be the same length as the number of chromosomes represented in the minimum distange genomic intervals (object)")
    }
    if(verbose)  pb <- txtProgressBar(min=0, max=length(indexList), style=3)
    segList <- vector("list", length(indexList))
    pkgs <- neededPkgs()
    j <- k <- NULL
    segList <- foreach(j=indexList, k=indexList2, featureData=fD, .packages=pkgs) %dopar%{
      narrowRangeForChromosome(object[j, ],
                               lrr.segs[k, ],
                               thr=thr,
                               verbose=FALSE,
                               fD=featureData)
    }
    if(verbose) close(pb)
    segs <- unlist(GRangesList(segList))
  } else {
    if(is(fD, "list")) {
      chroms <- paste("chr", sapply(fD, function(x) chromosome(x)[1]), sep="")
      fD <- fD[[match(as.character(chromosome(object))[1], chroms)]]
      chr <- as.character(chromosome(object)[1])
      if(paste("chr", chromosome(fD)[1], sep="") != chr) stop("The supplied GenomeAnnotatedDataFrame (fD) does not have the same chromosome as object")
    }
    segs <- narrowRangeForChromosome(object, lrr.segs, thr, verbose, fD=fD)
  }
  metadata(segs) <- metadata(lrr.segs)
  return(segs)
}


narrowRangeForChromosome <- function(md.range, cbs.segs, thr=0.9, verbose=TRUE, fD){
  md.range <- md.range[order(sampleNames(md.range), start(md.range)), ]
  mads <- pmax(values(md.range)$mindist.mad, .1)
  abs.thr <- abs(values(md.range)$seg.mean)/mads
  md.range2 <- md.range[abs.thr > thr, ]
  if(length(md.range2) < 1) return(md.range)

  cbs.segs <- cbs.segs[order(sampleNames(cbs.segs), start(cbs.segs)), ]
  o <- findOverlaps(md.range2, cbs.segs)
  j <- subjectHits(o)
  ## only consider the cbs segments that have an overlap
  if(!is.na(match("sample", colnames(cbs.segs)))) cbs.segs <- cbs.segs[, -match("sample", colnames(cbs.segs))]
  offspring.segs <- cbs.segs[j, ]
  sns <- unique(sampleNames(md.range2))
  chr <- chromosome(md.range)[1]
  rdlist <- list()
  for(j in seq_along(sns)){
    md <- md.range2[sampleNames(md.range2) == sns[j], ]
    of <- offspring.segs[sampleNames(offspring.segs)==sns[j], ]
    md.mad <- values(md)$mindist.mad
    md <- md[, match(colnames(values(of)), colnames(values(md)))]
    ## stack the ranges of the minimum distance segments and the offspring segments
    un <- unlist(GRangesList(list(md, of)))
    ## find the disjoint ranges
    disj <- disjoin(un)
    o <- findOverlaps(md, disj)
    ## which minimumdistance intervals are spanned by a disjoint interval
    r <- subjectHits(o)
    s <- queryHits(o)
    ## only keep the disjoint intervals for which a minimum distance segment is overlapping
    ##  (filters intervals that have a minimum distance of approx. zero)
    disj <- disj[r, ]
    elementMetadata(disj)$sample <- sampleNames(md)[s]
    elementMetadata(disj)$numberProbes <- 0L ## update later
    elementMetadata(disj)$seg.mean <- values(md)$seg.mean[s]
    elementMetadata(disj)$mindist.mad <- md.mad[s]
    rdlist[[j]] <- disj
  }
  rd <- unlist(GRangesList(rdlist))
  metadata(rd) <- metadata(md.range)
  frange <- makeFeatureGRanges(fD, metadata(rd)[["genome"]])
  cnt <- countOverlaps(rd, frange)
  elementMetadata(rd)$numberProbes <- cnt
  rd <- rd[numberProbes(rd) > 0L, ]
  mrd <- md.range[abs.thr <= thr, ]
  mrd <- mrd[, match(colnames(values(rd)), colnames(values(mrd)))]
  rd2 <- c(rd, mrd)
  rd2 <- rd2[order(sampleNames(rd2), start(rd2)), ]
  return(rd2)
}


##narrow <- function(object, lrr.segs, thr=0.9, mad.minimumdistance, verbose=TRUE){
##	if(!is(names(mad.minimumdistance), "character")) stop("mad.minimumdistance must be named")
##	##stopifnot(!is.null(names(mad.minimumdistance)))
##	ix <- match(sampleNames(object), names(mad.minimumdistance))
##	object$mindist.mad <- mad.minimumdistance[ix]
##	stopifnot("mindist.mad" %in% colnames(object))
##	lrr.segs <- lrr.segs[sampleNames(lrr.segs) %in% sampleNames(object), ]
##	if(length(unique(chromosome(object))) > 1){
##		if(verbose)
##			message("narrowing the ranges by chromosome")
##		indexList <- split(seq_len(nrow(object)), chromosome(object))
##		indexList2 <- split(seq_len(nrow(lrr.segs)), chromosome(lrr.segs))
##		stopifnot(all.equal(names(indexList), names(indexList2)))
##		if(verbose) {
##			pb <- txtProgressBar(min=0, max=length(indexList), style=3)
##		}
##		segList <- vector("list", length(indexList))
##		for(i in seq_along(indexList)){
##			if(verbose) setTxtProgressBar(pb, i)
##			j <- indexList[[i]]
##			k <- indexList2[[i]]
##			md.segs <- object[j, ]
##			lr.segs <- lrr.segs[k, ]
##			segList[[i]] <- narrowRangeForChromosome(md.segs, lr.segs, thr=thr, verbose=FALSE)
##			rm(md.segs, lr.segs); gc()
##		}
##		if(verbose) close(pb)
##		segs <- stack(RangedDataList(segList))
##		j <- match("sample", colnames(segs))
##		if(length(j) == 1) segs <- segs[, -j]
##	} else {
##		segs <- narrowRangeForChromosome(object, lrr.segs, thr, verbose)
##	}
##	rd.cbs <- RangedDataCBS(ranges=ranges(segs), values=values(segs))
##	return(rd.cbs)
##}
##
##
##narrowRangeForChromosome <- function(md.range, cbs.segs, thr=0.9, verbose=TRUE){
##	md.range <- md.range[order(sampleNames(md.range), start(md.range)), ]
##	cbs.segs <- cbs.segs[order(sampleNames(cbs.segs), start(cbs.segs)), ]
##	ir1 <- IRanges(start(md.range), end(md.range))
##	ir2 <- IRanges(start(cbs.segs), end(cbs.segs))
##	mm <- findOverlaps(ir1, ir2)
##	qhits <- queryHits(mm)
##	shits <- subjectHits(mm)
##	index <- which(sampleNames(md.range)[qhits] == sampleNames(cbs.segs)[shits])
##	if(length(index) > 0){
##		qhits <- qhits[index]
##		shits <- shits[index]
##	} else stop("no overlap")
##	##---------------------------------------------------------------------------
##	##
##	## only narrow the range if the minimum distance segment is
##	## bigger than some nominal value. Otherwise, we use the
##	## minimum distance range as is.
##	##
##	##
##	##---------------------------------------------------------------------------
##	##abs.thr <- abs(md.range$seg.mean)/md.range$mindist.mad > thr
##	mads <- pmax(md.range$mindist.mad, .1)
##	abs.thr <- abs(md.range$seg.mean)/mads > thr
##	## I1 is an indicator for whether to use the cbs start
##	deltaStart <- start(cbs.segs)[shits] > start(md.range)[qhits] & (start(cbs.segs)[shits] - start(md.range)[qhits] < 100e3)
##	deltaEnd <- end(cbs.segs)[shits] < end(md.range)[qhits]  & (end(md.range)[qhits] - end(cbs.segs)[shits] < 100e3)
##	I1 <- deltaStart & start(cbs.segs)[shits] >= start(md.range)[qhits] & start(cbs.segs)[shits] <= end(md.range)[qhits] & abs.thr[qhits]
##	## indicator of whether to use the cbs end
##	I2 <- deltaEnd & end(cbs.segs)[shits] <= end(md.range)[qhits] & end(cbs.segs)[shits] >= start(md.range)[qhits] & abs.thr[qhits]
##	st <- start(cbs.segs)[shits] * I1 + start(md.range)[qhits] * (1-I1)
##	en <- end(cbs.segs)[shits] * I2 + end(md.range)[qhits] * (1-I2)
##	st.index <- (cbs.segs$start.index[shits] * I1 + md.range$start.index[qhits]*(1-I1))
##	en.index <- (cbs.segs$end.index[shits] * I2 + md.range$end.index[qhits]*(1-I2))
##	## For each md.range range, there should only be one I1 that is TRUE
##	## If I1 and I2 are true, then a range is completely contained within the md.range segment
##	ids <- md.range$id[qhits]
##	##  |--------------|
##	## ---|--------|-----
##	## Becomes
##	##  |-|--------|---|
##	.i <- which(I1 & I2)
##	if(length(.i) > 0){
##		## new intervals
##		ir <- IRanges(st[.i], en[.i])
##		## remove intervals from ir1
##		## these intervals in ir1 must be bigger
##		ix <- subjectHits(findOverlaps(ir, ir1))
##		originalIntervals <- ir1[ix, ]
##		ids <- md.range$id[ix]
##		means <- md.range$seg.mean[ix]
##		chr <- chromosome(md.range)[ix]
##		mads <- md.range$mindist.mad[ix]
##		firstNew <- IRanges(start(originalIntervals), start(ir)-1)
##		endNew <- IRanges(end(ir)+1, end(originalIntervals))
##		ids <- rep(ids, each=3)
##		chr <- rep(chr, each=3)
##		means <- rep(means, each=3)
##		mads <- rep(mads, each=3)
##		res <- c(ir, firstNew, endNew)
##		## remove originalIntervals from the original set
##		ir1 <- ir1[-ix, ]
##		ids1 <- md.range$id[-ix]
##		chr1 <- chromosome(md.range)[-ix]
##		means1 <- md.range$seg.mean[-ix]
##		mads1 <- md.range$seg.mean[-ix]
##		irnew <- c(ir1, res)
##		ids2 <- c(ids1, ids)
##		chr2 <- c(chr1, chr)
##		mads2 <- c(mads1, mads)
##		means2 <- c(means1, means)
##		tmp <- RangedDataCBS(irnew,
##				     sampleId=ids2,
##				     chrom=chr2,
##				     seg.mean=means2,
##				     mindist.mad=mads2)
##	} else return(md.range)
##		##irnew <- irnew[order(start(irnew)), ]
####	index <- which(I1 & I2)-1
####	index <- index[index!=0]
####	index <- index[ids[index] == ids[index+1]]
#### 	ir3 <- IRanges(st[index], en[index])
####	ir4 <- IRanges(st[index+1], en[index+1])
####	##newRanges1 <- IRanges(st[index], st[index+1]-1)
####	##newRanges2 <- IRanges(st[index+1], en[index])
####	if(length(index) > 1){
####		originalEnd <- en[index]
####		originalStart <- st[index]
####		newEnd <- st[index+1]-1
####		newStart <- st[index+1]
####		en[index] <- newEnd
####		st <- c(st, newStart)
####		en <- c(en, originalEnd)
####		##en.index[index] <- st.index[index+1]-1
####	}
##	##split(I1, qhits)
##	##split(I2, qhits)
##	##stopifnot(sapply(split(st, qhits), function(x) all(diff(x) >= 0)))
####	nm <- apply(cbind(st.index, en.index), 1, function(x) length(x[1]:x[2]))
##	## keep segment means the same as the minimum distance
####	tmp <- RangedData(IRanges(st, en),
####			  ##id=sampleNames(md.range)[qhits],
####			  id=ids2,
####			  ##chrom=chromosome(md.range)[qhits],
####			  chrom=chr2,
####			  ##num.mark=nm,
####			  seg.mean=md.range$seg.mean[qhits],
####			  start.index=st.index,
####			  end.index=en.index,
####			  mindist.mad=md.range$mindist.mad[qhits])
##	##family=md.range$family[qhits])
##	##ranges.below.thr <- split(!abs.thr[qhits], qhits)
##	##ns <- sapply(ranges.below.thr, sum)
####	uid <- paste(tmp$id, start(tmp), tmp$chrom, sep="")
##	##duplicated(uid)
##	##stopifnot(!all(duplicated(uid)))
####	tmp <- tmp[!duplicated(uid), ]
##	## for each subject, the following must be true
####	index <- which(tmp$id[-nrow(tmp)] == tmp$id[-1])
##	##stopifnot(all(end(tmp)[index] < start(tmp)[index+1]))
##	res <- tmp[order(tmp$id, start(tmp)), ]
##	return(res)
##}

stackListByColIndex <- function(object, i, j){
	X <- vector("list", length(object))
	is.matrix <- is(object[[1]], "matrix") || is(object[[1]], "ff_matrix")
	if(is.matrix){
		for(k in seq_along(X)){
			X[[k]] <- object[[k]][, i, drop=FALSE]
		}
		X <- do.call("rbind", X)
	} else {
		is.array <- is(object[[1]], "array")
		stopifnot(length(j)==1)
		for(k in seq_along(X)){
			X[[k]] <- object[[k]][, i, j, drop=FALSE]
			dim(X[[k]]) <- c(nrow(X[[k]]), ncol(X[[k]])) ## drop 3rd dimension
		}
		X <- do.call("rbind", X)
	}
	return(X)
}


callDenovoSegments <- function(path="",
			       pedigreeData,
			       ext="",
			       featureData,
			       cdfname,
			       chromosome=1:22,
			       segmentParents,
			       mdThr=0.9,
			       prOutlierBAF=list(initial=1e-3, max=1e-1, maxROH=1e-3),
			       verbose=FALSE, genome=c("hg19", "hg18"), ...){
  pkgs <- c("GenomicRanges", "VanillaICE", "oligoClasses", "matrixStats", "MinimumDistance")
  genome <- match.arg(genome)
  if(!is(pedigreeData, "Pedigree")) stop("pedigreeData must be an object of class Pedigree")
  filenames <- file.path(path, paste(originalNames(allNames(pedigreeData)), ext, sep=""))
  ##obj <- read.bsfiles(filenames=filenames, path="", ext="")
  headers <- names(read.bsfiles(filenames[1], nrows=0))
  select <- match(c("SNP Name", "Allele1 - AB", "Allele2 - AB",
                    "Log R Ratio", "B Allele Freq"), headers)
##  ##labels <- setNames(c("Allele1", "Allele2", "LRR", "BAF"), keep[-1])
##  classes <- c(rep("character", 3), rep("numeric", 2))
##  header_info <- VanillaICE:::headerInfo(filenames[1], skip=10, sep=",",
##                                         keep=keep, labels=labels,
##                                         classes=classes)
##  obj <- read_beadstudio(filenames=filenames)
  obj <- read.bsfiles(filenames, select=select)
  if(missing(featureData)){
    trioSetList <- TrioSetList(lrr=integerMatrix(obj[, "lrr",], 100),
                               baf=integerMatrix(obj[, "baf",], 1000),
                               pedigreeData=pedigreeData,
                               chromosome=chromosome,
                               cdfname=cdfname,
                               genome=genome)
  } else {
    trioSetList <- TrioSetList(lrr=integerMatrix(obj[, "lrr",], 100),
                               baf=integerMatrix(obj[, "baf",], 1000),
                               pedigreeData=pedigreeData,
                               featureData=featureData,
                               chromosome=chromosome,
                               genome=genome)
  }
  isff <- is(lrr(trioSetList)[[1]], "ff")
  if(isff) pkgs <- c("ff", pkgs)
  md <- calculateMindist(lrr(trioSetList), verbose=verbose)
  mads.md <- mad2(md, byrow=FALSE)
  fns <- featureNames(trioSetList)
  md.segs <- segment2(object=md,
                      pos=position(trioSetList),
                      chrom=chromosome(trioSetList, as.list=TRUE),
                      verbose=verbose,
                      id=offspringNames(trioSetList),
                      featureNames=fns,
                      genome=genome,
                      ...)
  lrrs <- lrr(trioSetList)
  if(!segmentParents){
    ## when segmenting only the offspring,
    ## the trio names are the same as the sampleNames
    lrrs <- lapply(lrrs, function(x){
      dns <- dimnames(x)
      x <- x[, , 3, drop=FALSE]
      dim(x) <- c(nrow(x), ncol(x))
      dimnames(x) <- list(dns[[1]], dns[[2]])
      return(x)
    })
    id <- offspringNames(trioSetList)
  } else{
    id=trios(trioSetList)
  }
  pos <- position(trioSetList)
  lrr.segs <- segment2(object=lrrs,
                       pos=position(trioSetList),
                       chrom=chromosome(trioSetList, as.list=TRUE),
                       id=id, ## NULL if segmentParents is FALSE
                       verbose=verbose,
                       featureNames=fns, genome=genome, ...)
  md.segs2 <- narrowRanges(md.segs, lrr.segs, 0.9, mad.minimumdistance=mads.md, fD=Biobase::featureData(trioSetList))
  index <- splitIndicesByLength(seq_len(ncol(trioSetList)), 1)
  ##if(!getDoParRegistered()) registerDoSEQ()
  outdir <- ldPath()
  id <- sampleNames(trioSetList)
  object <- i <- NULL
  map.segs <- foreach(i=index,
                      .packages=pkgs, .combine="unlist") %dopar% {
                        MAP(object=trioSetList,
                            ranges=md.segs2,
                            id=id[i],
                            mdThr=mdThr,
                            prOutlierBAF=prOutlierBAF)
                      }
  return(map.segs)
}

make.unique2 <- function(names, sep="___DUP") make.unique(names, sep)
originalNames <- function(names){
	if(length(names) ==0) return(names)
	sep <- formals(make.unique2)[["sep"]]
	index <- grep(sep, names)
	if(length(index) > 0) names[index] <- sapply(names[index], function(x) strsplit(x, sep)[[1]][[1]])
	names
}

read.bsfiles2 <- function(path, filenames, sampleNames, z, marker.index,
			  lrrlist, baflist, featureNames){
	i <- seq_along(sampleNames)
	## this is simply to avoid having a large 'dat' object below.
	if(isPackageLoaded("ff")){
		NN <- min(length(sampleNames), 2)
		ilist <- splitIndicesByLength(i, NN)
		for(k in seq_along(ilist)){
			j <- ilist[[k]]
			sns <- sampleNames[j]
			dat <- read.bsfiles(path=path, filenames=filenames[j])
			dat <- dat[match(featureNames, rownames(dat)), , , drop=FALSE]
			l <- match(sns, colnames(baflist[[1]]))
			for(m in seq_along(marker.index)){
				M <- marker.index[[m]]
				baflist[[m]][, l, z] <- integerMatrix(as.matrix(dat[M, 2, ]), scale=1000)
				lrrlist[[m]][, l, z] <- integerMatrix(as.matrix(dat[M, 1, ]), scale=100)
			}
		}
		return(TRUE)
	} else {
		dat <- read.bsfiles(path=path, filenames=filenames)
	}
	return(dat)
}

stackRangedDataList <- function(...) {
	##object <- stack(RangedDataList(...))
	object <- GRangesList(list(...)[[1]])
	unlist(object)
	##j <- match("sample", colnames(object))
	##if(is.na(j))  object else object[, -j]
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~ The rest is old code that has been commented out.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##inferMissingRanges <- function(start, end){
##}
##
##inferNormalRangesFromPenn <- function(penn.object, vi.object){
##	chr <- chromosome(penn.object)
##	starts <- split(start(penn.object), chr)
##	ends <- split(end(penn.object), chr)
##	foreach(start=starts, end=ends) %do% inferMissingRanges(start=start,
##			      end=end)
##}


##shrinkTo <- function(x, x.0, DF.PRIOR){
##	DF <- ncol(x)-1
##	DF <- Ns-1
##	DF[DF < 1] <- 1
##	x.0 <- apply(x, 2, median, na.rm=TRUE)
##	x <- (x*DF + x.0*DF.PRIOR)/(DF.PRIOR + DF)
##	for(j in 1:ncol(x)) x[is.na(x[, j]), j] <- x.0[j]
##	return(x)
##}


##dna <- function(object) harmonizeDnaLabels(phenoData2(object[[1]])[, "DNA.Source", ])
##plate <- function(object) phenoData2(object[[1]])[, "Sample.Plate", ]

##readTable3 <- function(a=0.009){
##	## initialize with small value to avoid -Inf
##	results <- .C("calculateCHIT", a=a, M=array(0, dim=c(rep(5,6))))$M
##	## Make sure to transpose!
##	aperm(results)
##}

##arrangeSideBySide2 <- function(object1, object2){
##	grid.newpage()
##	lvp <- viewport(x=0,
##			y=0.05,
##			width=unit(0.50, "npc"),
##			height=unit(0.95, "npc"), just=c("left", "bottom"),
##			name="lvp")
##	pushViewport(lvp)
##	nfigs1 <- length(object1$condlevels[[1]])
##	nfigs2 <- length(object2$condlevels[[1]])
##	stopifnot(length(nfigs1) == length(nfigs2))
##	pushViewport(dataViewport(xscale=c(0,1), yscale=c(0.05,1), clip="on"))
##	object1$layout <- c(1, nfigs1)
##	print(object1, newpage=FALSE, prefix="plot1", more=TRUE)
##	upViewport(0)
##	lvp2 <- viewport(x=0.5,
##			 y=0.25,
##			 width=unit(0.50, "npc"),
##			 height=unit(0.95, "npc"), just=c("left", "bottom"),
##			 name="lvp2")
##	pushViewport(lvp2)
##	pushViewport(dataViewport(xscale=c(0,1), yscale=c(0.05,1), clip="on"))
##	object2$layout <- c(1, nfigs1)
##	print(object2, newpage=FALSE, prefix="plot2", more=TRUE)
##}


##read.bsfiles <- function(path="./", filenames, ext="", row.names=1,
##			 sep="\t",
##			 as.is=TRUE, header=TRUE,
##			 drop=FALSE, ...){
##	fnames <- file.path(path, paste(filenames, ext, sep=""))
##	stopifnot(all(file.exists(fnames)))
##	for(i in seq_along(filenames)){
##		cat(".")
##		tmp <- read.table(file.path(path, paste(filenames[i], ext, sep="")),
##				  row.names=row.names,
##				  sep=sep,
##				  header=header,
##				  as.is=as.is, ...)
##		if(i==1){
##			j <- grep("Log.R.Ratio", colnames(tmp))
##			k <- grep("B.Allele", colnames(tmp))
##			dat <- array(NA, dim=c(nrow(tmp), 2, length(filenames)))
##			if(!drop){
##				dimnames(dat) <- list(rownames(tmp),
##						      c("lrr", "baf"),
##						      basename(filenames))
##			}
##			##lrr.data <- matrix(NA, nrow(tmp), length(filenames))
##			##baf.data <- matrix(NA, nrow(tmp), length(filenames))
##		}
##		dat[, 1, i] <- tmp[, j]
##		dat[, 2, i] <- tmp[, k]
##	}
##	cat("\n")
##	return(dat)
##}

initializeLrrAndBafArrays <- function(dims, col.names, outdir, name=""){
	ldPath(outdir)
	if(name != ""){
		bafname <- paste(name, "baf", sep="_")
		lrrname <- paste(name, "lrr", sep="_")
	} else {
		bafname <- "baf"
		lrrname <- "lrr"
	}
	bafs <- initializeBigArray(bafname, dim=dims, vmode="integer")
	lrrs <- initializeBigArray(lrrname, dim=dims, vmode="integer")
	colnames(bafs) <- colnames(lrrs) <- col.names
	res <- list(baf=bafs, lrr=lrrs)
	return(res)
}

trioSetListExample <- function(){
	data(trioSetListExample)
	ad <- assayData(trioSetList)
	b <- lapply(ad[["BAF"]], integerArray, scale=1000)
	r <- lapply(ad[["logRRatio"]], integerArray, scale=100)
	ad2 <- AssayDataList(BAF=b, logRRatio=r)
	trioSetList@assayDataList <- ad2
	return(trioSetList)
}

neededPkgs <- function() c("oligoClasses", "Biobase", "MinimumDistance")

gcSubtractMatrix <- function(object, center=TRUE, gc, pos, smooth.gc=TRUE, ...){
	if(ncol(object) !=3) stop("Must pass one trio at a time.")
	cnhat <- matrix(NA, nrow(object), ncol(object))
	isna <- rowSums(is.na(object)) > 0
	if(any(isna)){
		i <- which(isna)
		gc <- gc[-i]
		pos <- pos[-i]
		if(smooth.gc)
			gc <- lowess(gc~pos, ...)$y
		object <- object[-i, , drop=FALSE]
	}
	X <- cbind(1, gc)
	## needs to be a big enough window such that we do not remove a true deletion/amplification
	index <- splitIndicesByLength2(x=seq_along(pos), lg=10000, MIN.LENGTH=2000)
	fitGCmodel <- function(X, Y){
		betahat <- solve(crossprod(X), crossprod(X, Y))
		yhat <- X %*% betahat
		resid <- Y-yhat
		resid
	}
	j <- NULL
	resid <- foreach(j=index, .combine="rbind") %do% fitGCmodel(X=X[j, ], Y=object[j, ])
	## ensure that the chromosome arm has the same median as in the original Y's
	if(center) resid <- resid+median(object,na.rm=TRUE)
	if(any(isna)) cnhat[-i, ] <- resid else cnhat <- resid
	cnhat
}


rescale2 <- function(x, l, u){
  y <- (x-min(x,na.rm=TRUE))/(max(x,na.rm=TRUE)-min(x,na.rm=TRUE))
  rescale(y, l, u)
}

dataFrameFromRange2 <- function(object, range, range.index, frame=0, nchar_sampleid=15){
	frange <- makeFeatureGRanges(featureData(object), genomeBuild(object))
	rm <- findOverlaps(range, frange, maxgap=frame)
	mm <- IRanges::as.matrix(rm)
	mm.df <- data.frame(mm)
	mm.df$featureNames <- featureNames(object)[mm.df$subject]
	marker.index <- mm.df$subject
	sample.index <- match(sampleNames(range), sampleNames(object))
	if(any(is.na(sample.index))) stop("sampleNames in RangedData do not match sampleNames in ", class(data), " object")
	sample.index <- unique(sample.index)
	obj <- object[marker.index, sample.index]
	mm.df$subject <- match(mm.df$featureNames, featureNames(obj))
	##
	## coersion to data.frame
	##
	df <- trioSet2data.frame(obj)
	chr <- unique(chromosome(object))
	oindex <- which(df$memberId=="offspring")
	findex <- which(df$memberId=="father")
	mindex <- which(df$memberId=="mother")
	memberid <- as.character(df$memberId)
	nchr <- min(nchar_sampleid, nchar(sampleNames(obj)[1]))
	oid <- substr(sampleNames(obj)[1], 1, nchr)
	fid <- substr(fatherNames(obj)[1], 1, nchr)
	mid <- substr(motherNames(obj)[1], 1, nchr)
	memberid[oindex] <- paste(oid, "(offspring)")
	memberid[findex] <- paste(fid, "(father)")
	memberid[mindex] <- paste(mid, "(mother)")
	memberid <- paste("chr", chr, ": ", memberid, sep="")
	memberid <- factor(memberid, levels=rev(unique(memberid)))
	df$memberId <- I(memberid)
	##df$range <- rep(i, nrow(df))##mm.df$query
	##dfList[[i]] <- df
	df$range <- range.index
	df
}

trioSet2data.frame <- function(from){
	cn <- lrr(from)[, 1, ]/100
	md.null <- is.null(mindist(from))
	if(!md.null) {
		md <- as.numeric(mindist(from))/100
		mdlabel <- "md"
	} else {
		mdlabel <- md <- NULL
	}
	labels <- c("father", "mother", "offspring", mdlabel)
	J <- length(labels)
	sns <- matrix(labels, nrow(cn), J, byrow=TRUE)
	sns <- as.character(sns)
	cn <- as.numeric(cn)
	y <- c(cn, md)
	bf <- as.numeric(baf(from)[, 1, ])/1000
	bf <- c(bf, rep(NA, length(md)))
	x <- rep(position(from)/1e6, J)
	is.snp <- rep(isSnp(from), J)
	df <- data.frame(x=x,
			 y=y,
			 baf=bf,
			 memberId=sns,
			 trioId=rep(sampleNames(from), length(y)),
			 is.snp=is.snp,
			 stringsAsFactors=FALSE)
	df$memberId <- factor(df$memberId, ordered=TRUE, levels=rev(labels))
	return(df)
}



segMeanAboveThr <- function(x, param){
  mads <- pmax(x$mindist.mad, minimum_MAD(param))
  abs(x$seg.mean/mads) > minimum_distance_threshold(param)
}

referenceIndex <- function(param) which(stateNames(param) == referenceState(param))
segmentLogLik <- function(log_emit){
  LLT <- apply(log_emit, c(2, 3), sum, na.rm=TRUE)
  ## copy number 2 prob is max(diploid not ROH, diploid ROH)
  LLT[, 3] <- pmax(LLT[, 3], LLT[, 4])
  ## remove diploid ROH state
  LLT <- LLT[, -4]
  LLT
}

loglik222 <- function(param, LLT){
  ## assume Pr(state_1,f | lambda) = Pr(state_2,m | lambda) = pi
  ## For offspring, we have Pr(state_1,o | state_1,f, state_1,m, DN=0, 1)
  ##    or 1/5 if DN=1
  ## if DN is 0 (not denovo), then many of the hidden
  ##  states should have  an epsilon probability of occurring.
  ##pi.offspring <- c(lookUpTable1(table1(penn), trio_state),  1/5)
  tau <- transitionProb(param)
  log.pi <- log(initialStateProb(param))
  state_index <- state(param)["222", ] + 1L
  lik <- diag(LLT[, state_index])

  prob_diploid <- log.pi[3]
  pr.nonmendelian <- prNonMendelian(param)
  pi.offspring <- c(table1(param)["222"], 1/5)
  lpr.offspring <- log(pi.offspring[1] * (1-pr.nonmendelian) + pi.offspring[2]*pr.nonmendelian)
  log.pi2 <- c(rep(prob_diploid, 2), lpr.offspring)
  sum(lik + log.pi2)
}

statesToEvaluate <- function(param, above_thr){
  nms <- stateNames(param)
  if(above_thr){
    x <- setNames(rep(TRUE, length(nms)), nms)
  } else {
    x <- setNames(rep(FALSE, length(nms)), nms)
    x[referenceIndex(param)] <- TRUE
  }
  x
}

.data_frame_loglik <- function(g){
  L <- length(g)
  data.frame(call=character(L),
             loglik=numeric(L),
             reference=numeric(L),
             LLR=numeric(L),
             row.names=names(g),
             stringsAsFactors=FALSE)
}

compute_loglik <- function(object, md_ranges, param){##, pr.nonmendelian,
  ##denovo.prev <- NULL
  g <- sort(md_ranges[numberProbes(md_ranges) > 2])
  uid <- unique(md_ranges$sample)
  if(length(uid) > 1) stop("md_ranges should contain only the offspring segments")
  hits <- findOverlaps(g, rowData(object))
  feature_index <- subjectHits(hits)
  results <- .data_frame_loglik(g)
  above_thr <- segMeanAboveThr(g, param)
  log_emit <- emissionArray(object, log_transform=TRUE,
                            epsilon=minimum_emission(param))
  trio_states <- state(param)
  state.prev <- NULL
  for(i in seq_along(g)){
    index <- feature_index[queryHits(hits) == i]
    LLT <- tryCatch(segmentLogLik(log_emit[index, , , drop=FALSE]),
                    error=function(e) NULL)
    if(is.null(LLT)) browser()
    evaluate <- statesToEvaluate(param, above_thr[i])
    states <- stateNames(param)[evaluate]
    result <- sapply(states, jointProb,
                     param=param,
                     state.prev=state.prev,
                     log.lik=LLT, USE.NAMES=FALSE)
    reference <- result[["222"]]
    loglik <- result[which.max(result)]
    results$loglik[i] <- round(loglik, 2)
    results$call[i] <- names(loglik)
    results$reference[i] <- reference
    state.prev <- trio_states[names(loglik), ]
  }
  results$call[!above_thr] <- NA
  results$LLR <- round(results$loglik-results$reference, 2)
  results
}

setSequenceLengths <- function(build, names){ ## names are unique(seqnames(object))
  sl <- getSequenceLengths(build)
  sl[match(unique(names), names(sl))]
}

.getArm <- function(chrom, pos, genome){
  if(is.integer(chrom)) chrom <- paste("chr", integer2chromosome(chrom), sep="")
  path.gap <- system.file("extdata", package="SNPchip")
  gaps <- readRDS(list.files(path.gap, pattern=paste("gap_", genome, ".rda", sep=""), full.names=TRUE))
  centromere.starts <- start(gaps)
  centromere.ends <- end(gaps)
  names(centromere.ends) <- names(centromere.starts) <- seqnames(gaps)
  centromere.starts <- centromere.starts[chrom]
  centromere.ends <- centromere.ends[chrom]
  chr.arm <- arm <- rep(NA, length(pos))
  arm[pos <= centromere.starts] <- "p"
  arm[pos >= centromere.ends] <- "q"
  ##arm <- ifelse(pos <= centromere.starts, "p", "q")
  chr.arm[!is.na(arm)] <- paste(chrom[!is.na(arm)], arm[!is.na(arm)], sep="")
  chr.arm
}


.checkOrder <- function(object, verbose=FALSE){
  d <- diff(order(chromosome(object), position(object)))
  if(any(d < 0)){
    if(verbose)
      warning("Object should be ordered by chromosome and physical position.\n",
              "Try \n",
              "> object <- order(object) \n")
    return(FALSE)
  }
  TRUE
}

isFF <- function(object){
  names <- ls(assayData(object))
  is(assayData(object)[[names[[1]]]], "ff") | is(assayData(object)[[names[[1]]]], "ffdf")
}


emissionArray <- function(object, log_transform=TRUE, epsilon=0.01){
  emitlist <- assays(object)
  if(log_transform){
    emitlist <- lapply(emitlist, function(x, epsilon) log(x+epsilon), epsilon=epsilon)
  }
  lemit_array <- array(NA, dim=c(nrow(object), 3, 6))
  lemit_array[, 1, ] <- emitlist[[1]]
  lemit_array[, 2, ] <- emitlist[[2]]
  lemit_array[, 3, ] <- emitlist[[3]]
  lemit_array
}
