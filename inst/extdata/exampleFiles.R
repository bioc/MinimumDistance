##---------------------------------------------------------------------------
##
## trioSetListExample: two chromosomes
##
##---------------------------------------------------------------------------
lrr_chr7and22 <- lrr(trioSetList)
baf_chr7and22 <- baf(trioSetList)
featureData_chr7and22 <- lapply(trioSetList, featureData)
save(lrr_chr7and22, file="~/Software/MinimumDistance/inst/extdata/lrr_chr7and22.rda")
save(baf_chr7and22, file="~/Software/MinimumDistance/inst/extdata/baf_chr7and22.rda")
save(featureData_chr7and22, file="~/Software/MinimumDistance/inst/extdata/featureData_chr7and22.rda")





##---------------------------------------------------------------------------
library(CleftExperimentData)
data(pedigreeCleft)
library(MinimumDistance)
pedigree1 <- pedigreeCleft[1,]
path <- "/thumper/ctsa/beaty/holger/txtfiles"
filenames <- file.path(path, paste(allNames(pedigree1), ".txt", sep=""))
file.copy(filenames[1], "~/Software/MinimumDistance/inst/extdata")

path <- "~/Software/MinimumDistance/inst/extdata"
fnames <- list.files(path, pattern=".txt", full.names=TRUE)
snames <- c("F", "M", "O")
cdfname <- "human610quadv1b"
for(i in seq_along(fnames)){
	dat <- read.table(fnames[i],
			  sep="\t",
			  header=TRUE,
			  as.is=TRUE)
	colnames(dat) <- c("Name", paste(snames[i], c("Log.R.Ratio", "B.Allele.Freq"), sep="."))

	if(i==1){
		featureAnnotation <- oligoClasses:::featureDataFrom(cdfname)
		fD <- featureAnnotation[order(featureAnnotation$chromosome, featureAnnotation$position), ]
		fD <- fD[order(fD$chromosome, fD$position), ]
		is.present <- sampleNames(fD) %in% dat$Name
		if(any(!is.present))
			fD <- fD[is.present, ]
	}
	index <- which(fD$chromosome==1)[1:1000]
	fns <- sampleNames(fD)[index]
	index <- match(fns, dat$Name)
	tmp <- dat[index, ]
	write.table(tmp,
		    file=fnames[i],
		    sep="\t",
		    row.names=FALSE,
		    quote=FALSE)
}


