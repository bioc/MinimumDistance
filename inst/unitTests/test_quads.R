test_RdsViews <- function(){
  library(data.table)
  library(VanillaICE)
  library(foreach)
  load_all("~/Software/Maher")
  pedigree_files <- file.path("~/Software/Maher/inst/extdata", c("quadsh.txt", "triosh.txt"))
  datdir <- "/dcs01/oncbio/rscharpf/maher"
  ##trace(PedigreeList, browser)
  pedlist <- do.call(c, sapply(pedigree_files, PedigreeList, path=datdir))
  ##trace(combineFamiliesInBothFiles, browser)
  pedlist2 <- combineFamiliesInBothFiles(pedlist)
  any(duplicated(unlist(pedlist))) ## some pedigrees have multiple generations

  views <- RdsViews(path=datdir, pedigree=pedlist2, cnvar="lrr", bafvar="baf")


  fgr <- readRDS(file.path(datdir, "feature_granges_hg19.rds"))
  assayList <- assays(views[1, ])
  ped <- pedigree(views[1,])
  me <- MinDistExperiment(assays=assayList,
                          rowData=fgr,
                          colData=setNames(DataFrame(ped), "filename"))
  ## check gender
  if(FALSE){
    ##index <- seq(1, nrow(me), 10)
    index <- which(chromosome(me)=="chrX")
    df <- data.frame(r=as.numeric(lrr(me)[index, ]),
                     b=as.numeric(baf(me)[index, ]),
                     pos=rep(start(me)[index], ncol(me)),
                     id=as.character(matrix(colnames(me), length(index), ncol(me), byrow=TRUE)))
    library(lattice)
    xyplot(r~pos | id, df, pch=20, cex=0.3, col="gray", layout=c(1,ncol(me)),
           panel=function(x,y,..., denovo_region){
             panel.xyplot(x, y, ...)
             panel.abline(h=0)
             panel.abline(h=-0.5, lty=2, col="gray")
           }, ylim=c(-2,2))
  }
  me <- subsetAndSort(me, seqlevels(me)[1:22])
  param <- MinDistParam()
  mdgr <- segment2(me, param=param)

  mdgr2 <- segment2(views[1, ], rowData=fgr, param=param)

  pedids <- pedigreeId(views)
  target_files <- file.path(path(views), paste0(pedids, "_mdgr.rds"))
  viewSapply(views[1:2], segment2, rowData=fgr, param=param,
             target_file=target_files[1:2])

  viewSapply(views[1:2], segment2, rowData=fgr, param=param,
             target_file=target_files[1:2])

  pedid <- gsub("_arraytv.rds", "", pedigree_id(mdgr))
  save(object, file=file.path(path(object), paste0(pedid, "_mdgr.rds")))

  target_file <-
  mdgr2 <- view(views[1, ],



  nMAD(param) <- 0.75
  mindist(mdgr) <- narrow2(mdgr, param)
  md_grl <- MAP2(me, mdgr, param)

  md_ranges <- md_grl[[1]]
  denovo <- md_ranges[isDenovo(md_ranges$call)]
  index <- subjectHits(findOverlaps(denovo, me, maxgap=20*width(denovo)))
  df <- data.frame(r=as.numeric(lrr(me)[index, ]),
                   b=as.numeric(baf(me)[index, ]),
                   pos=rep(start(me)[index], 3),
                   id=as.character(matrix(names(pedigree(me)), length(index), 3, byrow=TRUE)))
  library(lattice)
  xyplot(r~pos | id, df, pch=20, cex=0.3, col="gray", layout=c(1,3),
         panel=function(x,y,..., denovo_region){
           panel.xyplot(x, y, ...)
           panel.abline(v=c(start(denovo_region), end(denovo_region)))
         }, denovo_region=denovo)
  xyplot(b~pos | id, df, pch=20, cex=0.3, col="gray", layout=c(1,3),
         panel=function(x,y,..., denovo_region){
           panel.xyplot(x, y, ...)
           panel.abline(v=c(start(denovo_region), end(denovo_region)))
         }, denovo_region=denovo)
  xyplot(b~pos | id, df)

  ## do the entire pipeline
  MAP2(views[1, ], fgr, param)
}
