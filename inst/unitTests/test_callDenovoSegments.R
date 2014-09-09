test_callDenovoSegments <- function(){
  library(oligoClasses)
  foreach::registerDoSEQ()
  path <- system.file("extdata", package="MinimumDistance", mustWork=TRUE)
  fnames <- setNames(list.files(path, pattern=".txt", full.names=TRUE), c("father", "mother", "offspring1", "offspring2"))
  library(data.table)
  build <- "hg18"
  library("human610quadv1bCrlmm")



##
##  setClass("FileParam", representation(importfun="function", labels="character"))
##  ImportParam <- function(importfun=fread, labels=c("Name", "Log.R.Ratio", "B.Allele.Freq")){
##    new("FileParam", importfun=importfun, labels=labels)
##  }
##  setClass("AnnotationParam", representation(package="character", build="character"))
##  AnnotationParam <- function(package=character(), build="hg19"){
##    new("AnnotationParam", package=package, build=build)
##  }
##  setClass("FileViews", representation(path="character", pedigree="list"))
##  PedigreeList <- function(path="./", pedigree=list()){
##    new("FileViews", path=path, pedigree=pedigree)
##  }
##  FileViews <- function(import=ImportParam(), annotation=AnnotationParam(),
##                        pedigree=PedigreeList()){
##
##  }


##  views <- FileViews(annotation=AnnotationParam(build="hg18", package="human610quadv1bCrlmm"),
##                     pedigree=PedigreeList(path=path, pedigree=list(fnames)))
##
##  views <- FileViews(path=path, pedigree=list(fnames),
##                     fid="Name", cnvar="Log.R.Ratio", bafvar="B.Allele.Freq",
##                     headers=varLabels(),
##                     importfun=fread,
##                     annot_pkg=annot_pkg,
##                     build="hg18")
##  me <- MinDistExperiment(views)
  ##  map.segs <- callDenovoSegments(path=path,
##                                 ext="",
##                                 pedigreeData=ped,
##                                 cdfname="human610quadv1b",
##                                 chromosome=1,
##                                 segmentParents=FALSE,
##                                 genome="hg18")
##
}
