test_SnpArrayExperiment_pipeline <- function(){
  library(data.table)
  library(VanillaICE)
  library(Maher)
  datdir <- "/dcs01/oncbio/rscharpf/maher"
  atv.files <- list.files(datdir, full.names=TRUE, pattern="arraytv")
  ped <- fread("/dcs01/oncbio/rscharpf/maher/triosh.txt", header=FALSE)
  ##ped <- formatPedigree(ped)
  fgr2 <- readRDS(file.path(datdir, "feature_granges_hg19.rds"))

  views <- RdsViews(path=datdir, pedigree=ped, cnvar="lrr", bafvar="baf")
  assayList <- assays(views[1, ])
  me <- MinDistExperiment(assays=assayList,
                          rowData=fgr2,
                          pedigree=pedigree(views[1,]))
  me <- subsetAndSort(me, seqlevels(me)[1:22])
  param <- MinDistParam()
  path <- system.file("extdata", package="MinimumDistance")
  if(FALSE){
    md_granges <- segment2(me, param=DNAcopyParam())
    saveRDS(md_granges, file=file.path(path, "md_granges.rds"))
  } else md_granges <- readRDS(file.path(path, "md_granges.rds"))
  param <- PennParam()
  md_ranges <- MAP2(me, md_granges, param)
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
