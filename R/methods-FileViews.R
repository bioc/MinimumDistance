#' @export
FileViews <- function(path=character(), pedigree=list(),
                      cnvar=character(), bafvar=character(),
                      fid=character(),
                      importfun=function(){},
                      annot_pkg=character()){
  new("FileViews", path=path, pedigree=pedigree, cnvar=cnvar, bafvar=bafvar,
      fid=fid, importfun=importfun, annot_pkg=character())
}

.path <- function(object) object@path
filenames <- function(object) as.character(pedigree(object)[[1]])

setMethod("nrow", "FileViews", function(x) length(x@pedigree))
setMethod("length", "FileViews", function(x) length(x@pedigree))
setMethod("pedigree", "FileViews", function(object) object@pedigree)
cnvar <- function(object) object@cnvar
bafvar <- function(object) object@bafvar
fid <- function(object) object@fid

setMethod("show", "FileViews", function(object){
  cat("class 'FileViews'\n")
  cat("   No. pedigrees:", length(object), "\n")
  cat("   path to RDS files:", .path(object), "\n")
})

setMethod("[", "FileViews", function(x, i, j, ..., drop=FALSE){
  if(!missing(i)){
    x@pedigree <- pedigree(x)[i]
  }
  x
})

pednames <- function(object) names(pedigree(object)[[1]])

setMethod("assays", "FileViews", function(x, ...){
  if(nrow(x) > 1) stop("object contains several trios. Use [i, ] to subset the ith trio prior to calling assays")
  files <- file.path(.path(x), filenames(x))
  tmp <- lapply(files, x@importfun)
  tmp <- lapply(tmp, as.data.frame)
  cnindex <- grep(cnvar(x), colnames(tmp[[1]]))
  bindex <- grep(bafvar(x), colnames(tmp[[1]]))
  if(length(cnindex) != 1) stop("cnvar not correctly specified")
  if(length(bindex) != 1) stop("bafvar not correctly specified")
  lrrlist <- lapply(tmp, "[", cnindex)
  baflist <- lapply(tmp, "[", bindex)
  r <- setNames(do.call(cbind, lrrlist), pednames(x))
  b <- setNames(do.call(cbind, baflist), pednames(x))
  r <- as.matrix(r)
  b <- as.matrix(b)
  if(fid(x) != ""){
    findex <- grep(fid(x), colnames(tmp[[1]]))
    fid <- tmp[[1]][, findex]
    rownames(r) <- rownames(b) <- fid
  }
  snpArrayAssays(cn=r, baf=b)
})

##setGeneric("files", function(object) standardGeneric("files"))
##setMethod("FileViews", "FileViews", function(object){
##  file.path
##})


setMethod(MAP2, "FileViews", function(object, mdgr, param, rowData, ...){
  object <- object[1]
  assayList <- assays(object)
  ped <- pedigree(object)
  me <- MinDistExperiment(assays=assayList,
                          rowData=rowData,
                          colData=setNames(DataFrame(ped), "filename"))
  me <- subsetAndSort(me, seqlevels(me)[1:22])
  mdgr <- segment2(me, param=param)
  mindist(mdgr) <- narrow2(mdgr, param)
  mindist(mdgr) <- MAP2(me, mdgr, param)
})
