.test_SnpArrayExperiment_pipeline <- function(){
  ##
  ## Uses data not included in package
  ##
  library(data.table)
  library(VanillaICE)
  datdir <- system.file("extdata", package="VanillaICE")
  files <- list.files(datdir, pattern="FinalReport")
  pedlist <- list(ParentOffspring(id = "hapmap", father="12287_03",
                                  mother="12287_02",
                                  offspring="12287_01",
                                  filePaths=filePaths(views)[1:3]))
  ##
  ## Maher
  ##
  library(Maher)
  datdir <- "/dcs01/oncbio/rscharpf/maher"
  data(pedlist, package="Maher")
  atv.files <- list.files(datdir, full.names=TRUE, pattern="arraytv")
  ##ped <- fread("/dcs01/oncbio/rscharpf/maher/triosh.txt", header=FALSE)
  ##ped <- formatPedigree(ped)
  fgr2 <- readRDS(file.path(datdir, "feature_granges_hg19.rds"))
  views <- ArrayViews(rowRanges=fgr2,
                      filePaths=atv.files,
                      sample_ids=names(pedlist[[1]]))
  ##pedigree=ped, cnvar="lrr", bafvar="baf")
  ##assayList <- assays(views[1, ])
##  me <- MinDistExperiment(assays=assayList,
##                          rowRanges=fgr2,
  ##                          pedigree=pedigree(views[1,]))
  me <- MinDistExperiment(views, pedigree=pedlist[[1]])
  mdgr <- segment2(me, param=param)
  param <- MinDistParam()
  ##path <- system.file("extdata", package="MinimumDistance")
  ##param <- PennParam()
  md_ranges <- MAP2(me, mdgr, param)
  table(md_ranges$call)

  if(FALSE){
    denovo <- md_ranges[isDenovo(md_ranges$call)]
    index <- subjectHits(findOverlaps(denovo, me, maxgap=20*width(denovo)))
    ##
    ##xyplot(me, denovo)
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

    ##  "wget --timestamping \
    ## 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz'  \
    ## -O ~/hg19ToHg18.over.chain.gz"
    library(rtracklayer)
    system("gunzip ~/hg19ToHg18.over.chain.gz")
    chain <- import.chain("~/hg19ToHg18.over.chain")
    denovo.hg18 <- unlist(liftOver(denovo, chain))

    files <- list.files("~/Software/FilterVariants/inst/extdata", pattern="-06-03.rds$", full.names=TRUE)
    grl <- sapply(files, readRDS)
    findOverlaps(denovo.hg18, grl[[1]], maxgap=5e3)
    findOverlaps(denovo.hg18, grl[[2]], maxgap=5e3)
  }
}
