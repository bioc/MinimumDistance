quad_RdsViews <- function(){
  ##
  ## This example uses data not included in the package
  ##
  library(data.table)
  library(VanillaICE)
  library(foreach)
  pedigree_files <- file.path("~/Software/Maher/inst/extdata", c("quadsh.txt", "triosh.txt"))
  datdir <- "/dcs01/oncbio/rscharpf/maher"
  ##trace(PedigreeList, browser)
  pedlist <- do.call(c, sapply(pedigree_files, PedigreeList, path=datdir))
  pedlist2 <- combineFamiliesInBothFiles(pedlist)
  any(duplicated(unlist(pedlist))) ## some pedigrees have multiple generations
}
