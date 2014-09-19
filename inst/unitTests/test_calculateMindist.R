test_calculateMindist <- function(){
  library(oligoClasses)
  library(foreach)
  library(IRanges) ## stack method
  registerDoSEQ()
  data(trioSetListExample)
  mdlist <- calculateMindist(lrr(trioSetList))

  trioSet <- stack(trioSetList)
  md <- calculateMindist(lrr(trioSet))

  md1 <- do.call("rbind", mdlist)
  dimnames(md1) <- NULL
  dimnames(md) <- NULL
  checkTrue(identical(md, md1))

  ## array of log r ratios: F, M, O order
  ##   -- one marker, one trio:
  lrrArray <- array(NA, dim=c(1, 1, 3))
  lrrArray[, , ] <- c(0.2, 0.1, 0.12)
  checkEquals(as.numeric(calculateMindist(lrrArray)), 0.02, tolerance=0.0001)

  ## construct SnpArrayExperiment
  ## then, calculate min distance
  library(oligoClasses)
  path <- system.file("extdata", package="MinimumDistance")
  fnames <- list.files(path, pattern=".txt")
  ped <- Pedigree(data.frame(F=c("F.txt", "F.txt"),
                             M=c("M.txt", "M.txt"),
                             O=c("O.txt", "O1.txt")))




}
