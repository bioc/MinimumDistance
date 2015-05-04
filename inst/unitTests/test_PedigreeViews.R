##library(devtools)
##library(data.table)
##library(foreach)
##library(Maher)
##datdir <- "/dcs01/oncbio/rscharpf/maher"
##data(pedlist, package="Maher")
##pedlist <- lapply(pedlist, function(id, dir) file.path(dir, id), dir=datdir)
##fgr <- readRDS(file.path(datdir, "feature_granges_hg19.rds"))
##
####views <- ArrayViews(filePaths=unlist(pedlist),
####                    rowRanges=
##
##setClass("PedigreeView", contains="ArrayViews",
##         representation(pedigree='character'))
##
##
##PedigreeView <- function(..., pedigree=character()){
##  new("PedigreeView", ..., pedigree=pedigree)
##}
