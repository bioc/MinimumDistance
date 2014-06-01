test_callDenovoSegments <- function(){
  library(oligoClasses)
  foreach::registerDoSEQ()
  path <- system.file("extdata", package="MinimumDistance", mustWork=TRUE)
  fnames <- list.files(path, pattern=".txt")
  ped <- Pedigree(fatherIds=fnames[1], motherIds=fnames[2],
                  offspringIds=fnames[3])
  map.segs <- callDenovoSegments(path=path,
                                 ext="",
                                 pedigreeData=ped,
                                 cdfname="human610quadv1b",
                                 chromosome=1,
                                 segmentParents=FALSE,
                                 genome="hg18")

}
