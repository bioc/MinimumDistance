#' @export
MinDistGRanges <- function(mindist=GRangesList(), offspring=GRangesList(),
                           father=GRanges(),
                           mother=GRanges(),
                           pedigree_id=character(),
                           mad=numeric(), acf=numeric()){
  new("MinDistGRanges", mindist=mindist, offspring=offspring,
      father=father, mother=mother, pedigree_id=pedigree_id,
      mad=mad, acf=acf)
}


setMethod("names", "MinDistGRanges", function(x) x@id)
setMethod("mindist", "MinDistGRanges", function(object) object@mindist)

setReplaceMethod("mindist", c("MinDistGRanges", "GRangesList"), function(object, value) {
  object@mindist <- value
  object
})

setMethod("offspring", "MinDistGRanges", function(object) object@offspring)
setMethod("mother", "MinDistGRanges", function(object) object@mother)
setMethod("father", "MinDistGRanges", function(object) object@father)
setMethod("mad", "MinDistGRanges", function(x, center = median(x), constant = 1.4826, na.rm = FALSE, low = FALSE, high = FALSE) x@mad)
setMethod("acfs", "MinDistGRanges", function(x) x@acf)
pedigree_id <- function(object) object@pedigree_id

setMethod("show", "MinDistGRanges", function(object){
  cat("An object of class 'MinDistGRanges' \n")
  cat("  pedigree id:  ", pedigree_id(object), "\n")
  cat("  mindist:  GRangesList of length", length(offspring(object)), "\n")
  cat("  offspring:  GRangesList of length", length(offspring(object)), "\n")
  cat("  mother:  GRanges of length", length(mother(object)), "\n")
  cat("  father:  GRanges of length", length(father(object)), "\n")
  mads <- round(mad(object), 2)
  autocorrelations <- paste(round(acfs(object), 2), collapse=", ")
  cat("  MADs: ", mads, "\n")
  cat("  lag10 ACFs: ", autocorrelations, "\n")
})

#' @export
narrow2 <- function(object, param){
  .narrowGRangesList(object, param)
}

.narrowGRangesList <- function(object, param){
  offspring_grl <- offspring(object)
  mindist_grl <- mindist(object)
  mads <- mad(object)[names(mindist_grl)]
  mindist_grl2 <- foreach(md_gr=mindist_grl, offspr_gr=offspring_grl, md.mad=mads) %do%{
    .narrowMinDistGRanges(md_gr=md_gr, offspr_gr=offspr_gr, md.mad=md.mad, param=param)
  }
  setNames(GRangesList(mindist_grl2), names(mindist_grl))
}

.narrowMinDistGRanges <- function(md_gr, offspr_gr, md.mad, param){
  ## threshold is the number of mads
  keep <- segMeanAboveThr(mean=md_gr$seg.mean, mad=md.mad, nmad=nMAD(param))
  md.below.thr <- md_gr[!keep]
  md_gr <- md_gr[keep]
  if(length(md_gr) < 1) return(md_gr)
  offspr_gr <- subsetByOverlaps(offspr_gr, md_gr)
  disj <- disjoin(c(md_gr, offspr_gr))
  disj <- subsetByOverlaps(disj, md_gr)
  hits <- findOverlaps(md_gr, disj)
  ## which minimumdistance intervals are spanned by a disjoint interval
  j <- subjectHits(hits) ##r
  i <- queryHits(hits)   ##s
  ##disj$sample <- names(object)
  disj$sample <- md_gr$sample[1]
  disj$seg.mean <- NA
  disj$seg.mean[j] <- md_gr$seg.mean[i]
  disj$filename <- md_gr$filename[1]
  ##
  ## filtered minimum distance ranges (md < thr) will only be in the
  ## offspring segments
  ##
  mcols(md.below.thr) <- mcols(md.below.thr)[, -grep("numberProbes", colnames(mcols(md.below.thr)))]
  disj <- sort(c(disj, md.below.thr))
  return(disj)
}
