FilterParamMD <- function(state=trioStateNames(), seqnames=paste0("chr", 1:22), ...){
  param <- FilterParam(state=state, seqnames=seqnames, ...)
  as(param, "FilterParamMD")
}
