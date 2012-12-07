##setOldClass("ff_array")
##setOldClass("ff_matrix")
##setClass("LogRratioSet", contains="eSet")
##setClassUnion("matrixOrNULL", c("matrix", "NULL"))
##setClassUnion("arrayORff_array", c("array", "ff_array"))
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ Pedigree Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####setClass("Pedigree", representation(trios="data.frame", trioIndex="data.frame"))
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ TrioSet Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##setClass("TrioSet", contains="gSet",
##	 representation(fatherPhenoData="AnnotatedDataFrame",
##			motherPhenoData="AnnotatedDataFrame",
##			pedigree="Pedigree",
##			mindist="matrixOrNULL"))
##
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~ TrioSetList Class ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##setClass("TrioSetList", contains="gSetList",
##	 representation(pedigree="Pedigree",
##			##assayDataList="AssayData",
##			##phenoData="AnnotatedDataFrame",
##			fatherPhenoData="AnnotatedDataFrame",
##			motherPhenoData="AnnotatedDataFrame"))
##			##featureDataList="list",
##			##chromosome="integer"))


