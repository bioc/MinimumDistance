\name{xsegment}
\alias{cbs.segs}
\alias{md.segs}
\docType{data}
\title{Segmentation of log R ratios}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{
	data(cbs.segs)
	data(md.segs)
}
\format{
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
\dontrun{
	data("trioSetList")
	sampleNames(trioSetList) <- offspringNames(trioSetList)
	offspringNames(trioSetList) <- offspringNames(trioSetList)
	fatherNames(trioSetList) <- fatherNames(trioSetList)
	motherNames(trioSetList) <- motherNames(trioSetList)
	save(trioSetList, file="~/Projects/MinimumDistance/data/trioSetList.rda")
	## how the 'md.segs' object was created
	md.segs <- xsegment(object=trioSetList,
			    verbose=TRUE)
	save(md.segs,file="~/Projects/MinimumDistance/data/md.segs.rda")
	## how the 'cbs.segs' object was created
	cbs.segs <- xsegment(object=trioSetList,
			     segment.mindist=FALSE,
			     verbose=TRUE)
	save(cbs.segs,file="~/Projects/MinimumDistance/data/cbs.segs.rda")
}
data(cbs.segs)
}
\keyword{datasets}
