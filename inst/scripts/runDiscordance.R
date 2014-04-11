library(cacheSweave)
outdir <- "/local_data/r00/beaty"
setCacheDir(outdir)
library(snow)
library(oligoClasses)
##options(cluster=makeCluster(22, "SOCK"))
options(cluster=makeCluster(22, "SOCK"))
library(doSNOW)
registerDoSNOW(getCluster())
