generatorTransitionProbs <- function(chrom, position, build, ids, TAUP, tauMAX){
	S <- 6
	CHR <- paste("chr", integer2chromosome(chrom), sep="")
	chrarm <- .getArm(as.integer(chrom), position, build)
	chrarm <- factor(chrarm, unique(chrarm))
	sl <- getSequenceLengths(build)
	sl <- sl[unique(CHR)]

	toGRanges <- function(statePath, id){##, j){
		sp <- split(statePath, chrarm)
		rlist <- lapply(sp, Rle)
		pos <- split(position, chrarm)
		states <- x1 <- x2 <- x <- rl <- NULL
		starts <- foreach(rl=rlist, x=pos) %do% x[start(rl)]
		ends <- foreach(rl=rlist, x=pos) %do% x[end(rl)]
		statelist <- foreach(states=sp, rl=rlist) %do% states[start(rl)]
		chrlist <- split(CHR, chrarm)
		chrl <- foreach(chrom=chrlist, rl=rlist) %do% chrom[start(rl)]
		gr <- foreach(states=statelist,
			      x1=starts,
			      x2=ends,
			      rl=rlist,
			      chrom=chrl) %do% {
				      GRanges(chrom, IRanges(x1, x2),
					      numberProbes=width(rl),
					      sample=id,
					      state=states,
					      seqlengths=sl)
			      }
		gr <- unlist(GRangesList(gr))
	}
	if(missing(tauMAX)) tauMAX <- 1-5e-8
	d <- diff(position)
	p <- exp(-2*d/TAUP)
	minimum <- 1-1/((S-1)) + 0.01
	p <- pmax(p, minimum)
	p <- pmin(p, tauMAX)
	##
	## 1-(1-tau)*(S-1)*c1, c1=1.  can't be less than 0.8
	##
	if(any(d < 0)) p[d < 0] <- 0.8 ## different chromosomes
	tau <- p; rm(p)

	initialProb <- rep(1/S, S)
	##log.initial <- log(initialProb)
	tauC <- 1-tau ## probability that state t+1 != state t
	lP.N2N <- log(tau)
	lP.N2A <- log(tauC) ##probability normal -> altered
	lP.A2A <- lP.N2N
	lP.A2N <- lP.N2A ##probability altered -> normal
	lP.A2Astar <- lP.N2A ## probability altered -> different altered state
	if(length(CHR)==1){
		fr <- GRanges(rep(CHR, length(position)), IRanges(position, width=1))
	} else fr <- GRanges(CHR, IRanges(position, width=1))
	transitionProbs <- list(lP.N2N=lP.N2N,
				lP.N2A=lP.N2A,
				lP.A2A=lP.A2A,
				lP.A2N=lP.A2N,
				lP.A2Astar=lP.A2Astar,
				fr=fr)
	funs <- list(toGRanges=toGRanges,
		     tau=tau,
		     transitionProbs=transitionProbs)
	return(funs)
}

generatorMatrix <- function(rlist, blist, chr, center, snp.index, anyNP, ped){
	chr <- chr
	rlist <- rlist
	blist <- blist
	ped <- ped
	toMatrix <- function(j){
		rl <- lapply(rlist, function(x, j) x[, j, ,drop=TRUE], j=j)
		bl <- lapply(blist, function(x, j) x[, j, ,drop=TRUE], j=j)
		r <- do.call("rbind", rl)/100
		b <- do.call("rbind", bl)/1000
		ids <- as.character(trios(ped)[j, ])
		if(anyNP) b <- b[snp.index, , drop=FALSE]
		rownames(r) <- rownames(b) <- NULL
		colnames(r) <- colnames(b) <- ids
		if(center){
			if(all(as.integer(chr) <= 22)){
				meds <- colMedians(r, na.rm=TRUE)
				r <- sweep(r, 2, meds)
			} else{
				autosome.index <- which(chr <= 22)
				rauto <- r[autosome.index, ]
				meds <- colMedians(rauto, na.rm=TRUE)
				r[autosome.index, ] <- sweep(rauto, 2, meds)
			}
		}
		sapply(rlist, close)
		sapply(blist, close)
		list(r=r, b=b)
	}
	return(toMatrix)
}

generatorMatrix2 <- function(R, B, chr, center, snp.index, anyNP, ped){
	chr <- chr
	R <- R
	B <- B
	ped <- ped
	toMatrix <- function(j){
		open(R)
		open(B)
		r <- R[, j, , drop=TRUE]/100
		b <- B[, j, , drop=TRUE]/1000
##		rl <- lapply(rlist, function(x, j) x[, j, ,drop=TRUE], j=j)
##		bl <- lapply(blist, function(x, j) x[, j, ,drop=TRUE], j=j)
##		r <- do.call("rbind", rl)/100
##		b <- do.call("rbind", bl)/1000
		ids <- as.character(trios(ped)[j, ])
		if(anyNP) b <- b[snp.index, , drop=FALSE]
		rownames(r) <- rownames(b) <- NULL
		colnames(r) <- colnames(b) <- ids
		if(center){
			if(all(as.integer(chr) <= 22)){
				meds <- colMedians(r, na.rm=TRUE)
				r <- sweep(r, 2, meds)
			} else{
				autosome.index <- which(chr <= 22)
				rauto <- r[autosome.index, ]
				meds <- colMedians(rauto, na.rm=TRUE)
				r[autosome.index, ] <- sweep(rauto, 2, meds)
			}
		}
		close(R)
		close(B)
		list(r=r, b=b)
	}
	return(toMatrix)
}

generatorOverlapFeatures <- function(feature.granges){
	overlapFUN <- function(ranges){
		hits <- findOverlaps(ranges, feature.granges)
##		cnt <- countOverlaps(ranges, feature.granges)
##		mm <- findOverlaps(feature.granges, ranges)
		##return(mm=mm, cnt=cnt)
		return(hits)
	}
	return(overlapFUN)
}
