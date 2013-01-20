setOldClass("ff_array")
setOldClass("ff_matrix")
##setClass("LogRratioSet", contains="eSet")
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("arrayORff_array", c("array", "ff_array"))
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ Pedigree Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setClass("Pedigree", representation(trios="data.frame",
				    trioIndex="data.frame"))


##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ TrioSet Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setClass("TrioSet", contains="gSet",
	 representation(fatherPhenoData="AnnotatedDataFrame",
			motherPhenoData="AnnotatedDataFrame",
			pedigree="Pedigree",
			mindist="matrixOrNULL"))
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ TrioSetList Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setClass("TrioSetList", contains="gSetList",
	 representation(pedigree="Pedigree",
			##assayDataList="AssayData",
			##phenoData="AnnotatedDataFrame",
			fatherPhenoData="AnnotatedDataFrame",
			motherPhenoData="AnnotatedDataFrame"))
			##featureDataList="list",
			##chromosome="integer"))

setClass("Pedigree2", contains="DataFrame")
##Pedigree2 <- function(..., row.names=NULL, check.names=TRUE){
##	df <- DataFrame(..., row.names=row.names, check.names=check.names)
##	pdf <- as(df, "Pedigree2")
##}
##setClass("TrioSE", contains="SummarizedExperiment",
##	 representation(pedigree="Pedigree2"))#,
##	 prototype(
##		   assays=SimpleList(logRRatio=array(),
##		   BAF=array()),
##		   pedigree=Pedigree2()))



