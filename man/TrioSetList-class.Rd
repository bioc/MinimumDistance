\name{TrioSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{TrioSetList-class}
\alias{allNames,TrioSetList-method}
\alias{[,TrioSetList-method}
\alias{calculateMindist,TrioSetList-method}
\alias{computeBayesFactor,TrioSetList-method}
\alias{dims,TrioSetList-method}
\alias{fatherNames,TrioSetList-method}
\alias{initialize,TrioSetList-method}
\alias{lrr,TrioSetList-method}
\alias{mad.marker<-,TrioSetList,list-method}
\alias{mad.mindist<-,TrioSetList-method}
\alias{mad.sample<-,TrioSetList,matrix-method}
\alias{mindist<-,TrioSetList-method}
\alias{mindist,TrioSetList-method}
\alias{motherNames,TrioSetList-method}
\alias{ncol,TrioSetList-method}
\alias{nrow,TrioSetList-method}
\alias{offspringNames,TrioSetList-method}
\alias{pedigree,TrioSetList-method}
\alias{prune,TrioSetList,RangedDataCNV-method}
\alias{sampleSheet,TrioSetList-method}
\alias{stack,TrioSetList-method}
\alias{trios,TrioSetList-method}

\title{Class \code{"TrioSetList"}}
\description{
  A container storing pedigree information, as well as low-level
  statistical summaries used for copy number estimation: the log R
  ratios and B allele frequencies.  The list structure is organized by
  chromosome, where each element of the list is a \code{TrioSet} object.
}
\section{Objects from the Class}{

   A list of \code{\linkS4class{TrioSet}} objects.  Each element in the
 list is a \code{\linkS4class{TrioSet}} object for one
 chromosome. Currently, only autosomes are supported by this package, so
 instances of \code{\linkS4class{TrioSet}} typically have 22
 elements. See the function \code{\link{TrioSetList}} for instantiating
 an object of class \code{TrioSetList}.

}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} Each element is a
 \code{TrioSet} object. }
    \item{\code{pedigree}:}{Object of class \code{"Pedigree"}.  Contains
 information on the trio-relationships. }
    \item{\code{sampleSheet}:}{Object of class \code{"SampleSheet"}
 Contains meta-data for the experiment, such as chemistry plate. In
 future versions, this may be an \code{AnnotatedDataFrame} class.  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\linkS4class{AssayData}"}, by class "list", distance 2.
Class \code{"\linkS4class{list_or_ffdf}"}, by class "list", distance 2.
Class \code{"vectorORfactor"}, by class "list", distance 3.
}
\section{Methods}{
  \describe{
    \item{allNames}{\code{signature(object = "TrioSetList")}: Returns
 character vector of all the sample names.  Note that \code{sampleNames}
 of a \code{TrioSetList} object is not the same.  In particular,
 \code{sampleNames(object)} returns only the character vector of offspring
 ids which uniquely identify a trio.  Hence, a separate method,
 \code{allNames}, is supplied when all the sample ids in the data set
 are required.  Finally, note that the length of the vector returned by
 \code{allNames} is not necessarily a multiple of 3 as mothers and
 fathers with multiple offspring would be included in multiple trios.}
    \item{[}{\code{signature(x = "TrioSetList")}: The "[" method is
 applied via \code{lapply} to each \code{TrioSet} element. }}
    \item{calculateMindist}{\code{signature(object = "TrioSetList")}:
    See \code{\link{calculateMindist}}.}
  \item{computeBayesFactor}{\code{signature(object = "TrioSetList")}:
    See \code{\link{computeBayesFactor}}.}
  \item{dims}{\code{signature(object = "TrioSetList")}: Return
 dimensions of the assayData elements}
    \item{fatherNames}{\code{signature(object = "TrioSetList")}:
 character vector of father ids.  See \code{\link{fatherNames}} for usage.}
    \item{initialize}{\code{signature(.Object = "TrioSetList")}: see
      also \code{\link{TrioSetList}}}
    \item{lrr}{\code{signature(object="TrioSetList")} Returns list of
 log R ratios for each chromosome.}
    \item{mad.marker<-}{\code{signature(x = "TrioSetList", value =
 "matrix")}: See also \code{\link{mad2}}}
    \item{mad.mindist<-}{\code{signature(x = "TrioSetList")}: See also \code{\link{mad2}}}
    \item{mad.sample<-}{\code{signature(x = "TrioSetList", value = "list")}: See also \code{\link{mad2}} }
    \item{mindist<-}{\code{signature(object = "TrioSetList")}: Adds the
 minimum distance matrix to each \code{TrioSet} element in the
 \code{TrioSetList} object. See \code{\link{mindist<-}} for usage.}
    \item{mindist}{\code{signature(object = "TrioSetList")}: See also \code{link{mindist}}.}
%    \item{minimumDistance}{\code{signature(object = "TrioSetList")}:}
    \item{motherNames}{\code{signature(object = "TrioSetList")}: See
 also \code{\link{motherNames}} for usage.}
    \item{ncol}{\code{signature(x = "TrioSetList")}: the number of
 trios, or equivalently the number of offspring.}
    \item{nrow}{\code{signature(x = "TrioSetList")}: The number of features.}
    \item{offspringNames}{\code{signature(object = "TrioSetList")}: See
 also \code{\link{offspringNames}} for usage.  Note that the offspring
 ids uniquely identify a trio, and the method \code{sampleNames} will
 return the same result.}
    \item{pedigree}{\code{signature(object = "TrioSetList")}: Accessor
 for pedigree information.  See also \code{\linkS4class{Pedigree}}}
    \item{prune}{\code{signature(object = "TrioSetList", ranges =
 "RangedDataCNV")}: See \code{\link{prune}} for usage.}
    \item{sampleSheet}{\code{signature(object = "TrioSetList")}:
 Accessor for metadata of the experiment. See \code{\linkS4class{SampleSheet}}.}
    \item{stack}{\code{signature(object = "TrioSetList")}: creates a
 \code{TrioSet} object from a \code{TrioSetList} object. See
 \code{\link{stack}} for usage.}
    \item{trios}{\code{signature(object = "TrioSetList")}: Returns a
 \code{data.frame} of the trios. Each row in the {data.frame} contains
 the sample identifiers for the father, mother, and offspring.}
}

\author{
R. Scharpf
}
\seealso{
  \code{\linkS4class{TrioSet}}, \code{\linkS4class{Pedigree}}

}
\examples{
showClass("TrioSetList")
data(trioSetListExample)
pedigree(trioSetList)
sampleSheet(trioSetList)
trioSet <- stack(trioSetList)
}
\keyword{classes}
