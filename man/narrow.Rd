\name{narrow}
\alias{narrow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Narrow the minimum distance ranges by the segmentation of the offspring
  copy number estimates}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
narrow(md.range, cbs.segs, thr, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{md.range}{
%%     ~~Describe \code{md.range} here~~
}
  \item{cbs.segs}{
%%     ~~Describe \code{cbs.segs} here~~
}
  \item{thr}{
%%     ~~Describe \code{thr} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (md.range, cbs.segs, thr, verbose = TRUE)
{
    i <- which(sampleNames(cbs.segs) \%in\% sampleNames(md.range))
    cbs.segs <- cbs.segs[i, ]
    i <- which(chromosome(cbs.segs) \%in\% chromosome(md.range))
    if (length(i) > 0 & length(i) < nrow(cbs.segs)) {
        cbs.segs <- cbs.segs[i, ]
    }
    chroms <- unique(chromosome(md.range))
    rdN <- list()
    if (verbose) {
        message("Narrowing the ranges with segment mean above ",
            thr, " if the offspring has a start/stop in the region")
        pb <- txtProgressBar(min = 0, max = length(chroms), style = 3)
    }
    for (i in seq_along(chroms)) {
        if (verbose)
            setTxtProgressBar(pb, i)
        j <- chroms[i]
        md <- md.range[chromosome(md.range) == j, ]
        md <- md[order(sampleNames(md), start(md)), ]
        jj <- which(chromosome(cbs.segs) == j)
        cbs <- cbs.segs[jj, ]
        cbs <- cbs[order(sampleNames(cbs), start(cbs)), ]
        ir1 <- IRanges(start(md), end(md))
        ir2 <- IRanges(start(cbs), end(cbs))
        mm <- findOverlaps(ir1, ir2)
        qhits <- queryHits(mm)
        shits <- subjectHits(mm)
        index <- which(sampleNames(md)[qhits] == sampleNames(cbs)[shits])
        if (length(index) > 0) {
            qhits <- qhits[index]
            shits <- shits[index]
        }
        else stop("no overlap")
        abs.thr <- abs(md$seg.mean) > thr
        I1 <- start(cbs)[shits] >= start(md)[qhits] & start(cbs)[shits] <=
            end(md)[qhits] & abs.thr[qhits]
        I2 <- end(cbs)[shits] <= end(md)[qhits] & end(cbs)[shits] >=
            start(md)[qhits] & abs.thr[qhits]
        st <- start(cbs)[shits] * I1 + start(md)[qhits] * (1 -
            I1)
        en <- end(cbs)[shits] * I2 + end(md)[qhits] * (1 - I2)
        st.index <- (cbs$start.index[shits] * I1 + md$start.index[qhits] *
            (1 - I1))
        en.index <- (cbs$end.index[shits] * I2 + md$end.index[qhits] *
            (1 - I2))
        ids <- md$id[qhits]
        index <- which(I1 & I2) - 1
        index <- index[index != 0]
        index <- index[ids[index] == ids[index + 1]]
        if (length(index) > 1) {
            en[index] <- st[index + 1] - 1
            en.index[index] <- st.index[index + 1] - 1
        }
        stopifnot(sapply(split(st, qhits), function(x) all(diff(x) >=
            0)))
        nm <- apply(cbind(st.index, en.index), 1, function(x) length(x[1]:x[2]))
        tmp <- RangedData(IRanges(st, en), id = sampleNames(md)[qhits],
            chrom = chromosome(md)[qhits], num.mark = nm, seg.mean = md$seg.mean[qhits],
            start.index = st.index, end.index = en.index, mindist.mad = md$mindist.mad[qhits])
        uid <- paste(tmp$id, start(tmp), tmp$chrom, sep = "")
        stopifnot(!all(duplicated(uid)))
        tmp <- tmp[!duplicated(uid), ]
        index <- which(tmp$id[-nrow(tmp)] == tmp$id[-1])
        stopifnot(all(end(tmp)[index] < start(tmp)[index + 1]))
        rdN[[i]] <- tmp[order(tmp$id, start(tmp)), ]
    }
    if (verbose)
        close(pb)
    if (length(rdN) > 1) {
        rdL <- stack(RangedDataList(rdN))
        rd <- rdL[, -grep("sample", colnames(rdL))]
        rdCbs <- RangedDataCBS(ranges = ranges(rd), values = values(rd))
    }
    else rdCbs <- RangedDataCBS(ranges = ranges(rdN[[1]]), values = values(rdN[[1]]))
    return(rdCbs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
