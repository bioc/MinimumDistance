\name{RangedData-classes}
\Rdversion{1.1}
\docType{class}
\alias{RangedDataCNV-class}
\alias{RangedDataCBS-class}
\alias{RangedDataHMM-class}
\alias{chromosome,RangedDataCNV-method}
\alias{coverage2,RangedDataCNV-method}
\alias{sampleNames,RangedDataCNV-method}
\alias{state,RangedDataCNV-method}
\alias{todf,RangedDataCNV,ANY-method}
\alias{trioNames,RangedDataCNV-method}

\title{Classes in MinimumDistance for data on ranges}
\description{

  \code{RangedDataCNV} is a class extending the virtual class
  \code{RangedDataCopyNumber}. \code{RangedDataCopyNumber} extends
  \code{\linkS4class{RangedData}}.

  \code{RangedDataCBS} extends \code{\linkS4class{RangedDataCNV}} and is
  useful for storing ranges from segmentation algorithms such as
  circular binary segmentation.  In particular, the columns 'chrom',
  'id', and 'num.mark' are required for instances of the class.

  \code{RangedDataHMM} extends \code{\linkS4class{RangedDataHMM}} and is
  useful for storing ranges from hidden Markov models such as PennCNV or
  VanillaICE.  A column labeled 'state' is required for the class.

}
\section{Objects from the Class}{

See \code{RangedDataCBS} and \code{RangedDataHMM} constructors.

}
\section{Slots}{
  \describe{
    \item{\code{ranges}:}{Object of class \code{"RangesList"} ~~ }
    \item{\code{values}:}{Object of class \code{"SplitDataFrameList"} ~~ }
    \item{\code{elementType}:}{Object of class \code{"character"} ~~ }
    \item{\code{elementMetadata}:}{Object of class \code{"DataTableORNULL"} ~~ }
    \item{\code{metadata}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{RangedDataCNV}"}, directly.
%Class \code{"\linkS4class{RangedDataCopyNumber}"}, by class
%"RangedDataCNV", distance 2.
Class \code{RangedDataCopyNumber}, by class "RangedDataCNV", distance 2.
Class \code{"\linkS4class{RangedData}"}, by class "RangedDataCNV", distance 3.
Class \code{"\linkS4class{DataTable}"}, by class "RangedDataCNV", distance 4.
Class \code{"\linkS4class{List}"}, by class "RangedDataCNV", distance 4.
Class \code{"\linkS4class{DataTableORNULL}"}, by class "RangedDataCNV", distance 5.
Class \code{"\linkS4class{Vector}"}, by class "RangedDataCNV", distance 5.
Class \code{"\linkS4class{Annotated}"}, by class "RangedDataCNV", distance 6.
}

\section{Methods}{
  \describe{
    \item{coverage2}{\code{signature(object = "RangedDataCNV")}: ... }
    \item{prune}{\code{signature(object = "TrioSetList", ranges = "RangedDataCNV")}: ... }
    \item{prune}{\code{signature(object = "TrioSet", ranges = "RangedDataCNV")}: ... }
    \item{state}{\code{signature(object = "RangedDataCNV")}: ... }
    \item{todf}{\code{signature(object = "RangedDataCNV", range = "ANY")}: ... }
	 }
}

\author{
R. Scharpf
}

\seealso{

  See \code{\linkS4class{RangedData}} for a detailed description and
  additional methods.

}
\examples{showClass("RangedDataCNV")}
\keyword{classes}
