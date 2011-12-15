\name{TrioSet-class}
\Rdversion{1.1}
\docType{class}
\alias{TrioSet-class}
\alias{[,TrioSet-method}
\alias{baf<-,TrioSet-method}
\alias{baf<-}
\alias{calculateMindist,TrioSet-method}
\alias{dim,TrioSet-method}
\alias{fatherNames<-,TrioSet,character-method}
\alias{fatherNames<-}
\alias{fatherNames,TrioSet-method}
\alias{initialize,TrioSet-method}
\alias{lrr<-}
\alias{mad.marker<-,TrioSet,matrix-method}
\alias{mad.marker<-,TrioSet,numeric-method}
\alias{mad.marker<-,TrioSet,NULL-method}
\alias{mad.marker}
\alias{mad.marker<-}
\alias{mad.mindist<-,TrioSet,numeric-method}
\alias{mad.mindist<-}
\alias{mad.sample<-,TrioSet,matrix-method}
\alias{mad.sample<-}
\alias{mindist<-,TrioSet,matrix-method}
\alias{mindist<-,TrioSet,ff_matrix-method}
\alias{mindist<-,TrioSet,NULL-method}
\alias{mindist,TrioSet-method}
\alias{mindist}
\alias{mindist<-}
\alias{motherNames<-,TrioSet,character-method}
\alias{motherNames,TrioSet-method}
\alias{motherNames<-}
\alias{offspringNames<-,TrioSet,character-method}
\alias{offspringNames,TrioSet-method}
\alias{offspringNames<-}
\alias{order,TrioSet-method}
\alias{phenoData2,TrioSet-method}
\alias{phenoData2}
\alias{prune,TrioSet,RangedDataCNV-method}
\alias{updateObject,TrioSet-method}
\alias{xyplot,formula,TrioSet-method}

\title{Class \code{"TrioSet"}}
\description{

    A TrioSet is a container for storing high throughput assay data and
  metadata from genotyping arrays when the study design is case-parent
  trios. In our application, de novo copy number alterations in affected
  offspring were of the primary interest.  Examination of the joint
  distribution of the log R ratios and B allele frequencies across
  members in the trios was a motivating factor in developing the TrioSet
  class. In particular, the assay data elements for a \code{TrioSet}
  object are 3-dimensional arrays: marker x trio x individual.
  Typically, a TrioSet instance stores the data for a single chromosome
  and an instance of TrioSetList is a list of TrioSets.  While having a
  single TrioSet for the entire dataset would simplify the classes of
  objects defined in this package, multiple arrays with thousands of
  trios and roughly a million markers are impractical on most machines.

}

%\section{Slots}{
%  \describe{
%    \item{\code{phenoData2}:}{Object of class \code{"array"} ~~ }
%    \item{\code{mindist}:}{Object of class \code{"matrix"} ~~ }
%    \item{\code{mad}:}{Object of class \code{"matrix"} ~~ }
%    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
%    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
%    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
%    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} ~~ }
%    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
%    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
%    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
%  }
%}
%\section{Extends}{
%Class \code{"LogRratioSet"}, directly.
%Class \code{"\linkS4class{eSet}"}, by class "LogRratioSet", distance 2.
%Class \code{"\linkS4class{VersionedBiobase}"}, by class "LogRratioSet", distance 3.
%Class \code{"\linkS4class{Versioned}"}, by class "LogRratioSet", distance 4.
%}

\section{Objects from the Class}{

  Objects can be created by calls of the form :

  \code{new("TrioSet", logRRatio, BAF, phenoArray, mindist, mad, ...)}.

}

\section{Subsetting}{
  \describe{

  As the \code{TrioSet} class is an extension of \code{eSet}, the
  subsetting is similar. One important difference is that the assay date
  elements are 3 dimensional arrays.  While \code{k} is not a formal
  argument in the generic for \code{"["}, \code{k} can be passed to the
  \code{"["} method for \code{TrioSet} objects for subsetting the 3rd
  dimension of the assay data.

  \item{}{\code{object[i, j, drop]}:
    \code{i} selects features, \code{j} selects trios, and \code{k}
    (though not part of the generic) selects for the individual in a
    trio.  Valid values for \code{k} are 1 (selects father), 2 (selects
    mother), or 3 (selects offspring).
  }

  }

}

\section{Acessors}{
  \describe{
    The \code{object} in the accessor descriptions that follow is a
    \code{TrioSet}:

    \item{}{\code{baf(object) <- value}: assign B allele frequencies.
      \code{value} is a 3-dimensional array (feature x trio x sample).}

    \item{}{\code{dim(object)}: Returns the dimension of
    the assay data elements. Each assay data element is a three
    dimensional array with dimensions for features, trios, and sample,
    respectively.}

    \item{}{\code{fatherNames(object) <- value}: Assign character string
    of father sample names.}

  \item{}{\code{fatherNames(object)}: Return character string of father
    identifiers.}

  \item{}{\code{lrr(object) <- value}: assign log R ratios to assayData.
      \code{value} is a 3-dimensional array (feature x trio x sample).}

  \item{}{\code{mad.marker(object) <- value}: Assign robust estimates of
  the log R ratio standard deviation for each marker (estimated across
  samples) to the \code{featureData}. \code{value} can be numeric or a
  matrix.}

\item{}{\code{mad.mindist(object) <- value}: Assign robust estimates of
  the minimum distance standard deviation for each trio (estimated across
  all markers) to the \code{phenoData}. \code{value} is a numeric vector.}


\item{}{\code{mad.sample(object) <- value}: Assign robust estimates of
  the log R ratio standard deviation for each sample (estimated across
   markers). \code{value} is a numeric vector.}

 \item{}{\code{mindist(object) <- value}: \code{value} is the minimum
   distance matrix. See \code{\link{calculateMindist}}}

 \item{}{\code{mindist(object)}: Accessor for the minimum distance matrix.}

 \item{}{\code{motherNames(object) <- value}: Assign a character string
 for the mother identifiers }

   \item{}{\code{motherNames(object)} Return a character vector of
     mother identifiers}

   \item{}{\code{offspringNames(object) <- value}: Assign a character
     vector of offspring identifiers}

   \item{}{\code{offspringNames(object)}:  Retrieve character vector of
 offspring identifiers.  Note that the result will identical to
 \code{sampleNames(object)} as the offspring identifiers uniquely
 identify a trio}

   \item{}{\code{order(object)}: order \code{TrioSet} object by
   chromosome and physical position}

    \item{}{\code{phenoData2(object)}: metadata for the samples stored
      as an array.  'sampleNames' is a required column header.}

  }

}

\section{Computing the minimum distance}{
  \describe{

    The \code{TrioSet} container stores the log R ratios for each member
    in the trio.  The 'minimum distance' is the minimum signed absolute
    difference of the parental log R ratios and the offspring log R
    ratios. Specifically, let |O-F| denote the absolute difference in
    the log R ratios comparing offspring to father and |O-M| the
    absolute difference in the log R ratios comparing offspring to
    mother.  The minimum distance at a marker is the signed minimum of
    |O-M| and |O-F|.  Plotted for all markers, deflections from zero can
    indicate a difference in the log R ratio of the offspring and either
    parent (possibly de novo).  For example, a positive minimum distance
    suggests that the log R ratio from the offspring is greater than the
    log R ratio of either parent.  The minimum distance can be computed
    from a \code{object} of class \code{TrioSet} by:
%
\item{}{\code{calculateMindist(object, ..., verbose=TRUE)}:

  The returned value is a \code{list}. Each element in the list contains
  a matrix of the minimum distances from one chromosome (dimension:
  number features x number trios).  Note that the \code{list} can be assigned
  back to the \code{TrioSetList} object using the replacement method
  \code{mindist<-}.  See the example.}

}

}

\section{Visualization}{
  \describe{
    \item{}{\code{xyplot(x, data, ...)}       }

  }
}

\author{R. Scharpf}


\seealso{
	\code{\linkS4class{TrioSetList}}
      }

\examples{
	showClass("TrioSet")
	## instantiate a TrioSet with an array of log R ratios (logRR) and B
	## allele frequencies
	logRR <- array(rnorm(10*2*3), dim=c(10, 2, 3))
	dimnames(logRR) <- list(paste("marker", 1:10, sep=""), paste("trio", 1:2, sep=""), c("F", "M", "O"))
	BAFs <- array(runif(10*2*3, 0, 1), dim=dim(logRR), dimnames=dimnames(logRR))
	trioSet <- new("TrioSet", logRRatio=logRR, BAF=BAFs)
	## log R ratios for first five markers for all trios
	lrr(trioSet)[1:5, , ]
	## log R ratios for first five markers and the first trio
	lrr(trioSet)[1:5, 1, ]
	## log R ratios for only the offspring in the first trio
	lrr(trioSet)[1:5, 1, "O"]
	## similarly for B allele frequencies
	baf(trioSet)[1:5, 1, "O"]
	if(require("ff")){
		## for large datasets, one can generate ff arrays for the assay data
		## rather than store these in memory
		##
		## First, indicate the directory to store the ff objects
		\dontrun{
			ldPath(".")
		}
		logRR.ff <- oligoClasses:::initializeBigArray(name="lrr", dim=c(10,2,3), vmode="double")
		## after initializion, one can read subsets of the data and write
		## to disk by the "[<-" assignment method
		## For example, here we write the data for 10 markers to disk
		logRR.ff[1:10, , ] <- logRR[1:10, , ]
		## Similarly for BAFs
		BAF.ff <- oligoClasses:::initializeBigArray(name="baf", dim=c(10,2,3), vmode="double")
		BAF.ff[1:10, , ] <- BAFs[1:10, , ]
		trioSet.ff <- new("TrioSet", logRRatio=logRR.ff, BAF=BAF.ff)
	}
}
\keyword{classes}
