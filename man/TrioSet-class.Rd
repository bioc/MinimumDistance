\name{TrioSet-class}
\Rdversion{1.1}
\docType{class}
\alias{TrioSet-class}
\alias{[,TrioSet-method}
\alias{[,ff_array-method}
\alias{allNames,TrioSet-method}
\alias{baf,TrioSet-method}
\alias{baf<-,TrioSet-method}
\alias{baf<-,TrioSet,array-method}
\alias{baf<-,TrioSet,ff_array-method}
\alias{baf<-}
\alias{calculateMindist,TrioSet-method}
\alias{coerce,TrioSet,TrioSetList-method}
\alias{coerce,TrioSet,data.frame-method}
\alias{computeBayesFactor,TrioSet-method}
\alias{dim,TrioSet-method}
\alias{fatherNames<-,TrioSet,character-method}
\alias{fatherNames<-}
\alias{fatherNames,TrioSet-method}
\alias{fatherPhenoData,TrioSet-method}
\alias{motherPhenoData,TrioSet-method}
\alias{offspringPhenoData}
\alias{offspringPhenoData,TrioSet-method}
\alias{initialize,TrioSet-method}
\alias{lrr<-}
\alias{lrr,TrioSet-method}
\alias{lrr<-,TrioSet-method}
\alias{mindist<-,TrioSet,matrix-method}
\alias{mindist,TrioSet-method}
\alias{mindist}
\alias{mindist<-}
\alias{motherNames,TrioSet-method}
\alias{ncol,TrioSet-method}
\alias{offspringNames,TrioSet-method}
\alias{order,TrioSet-method}
\alias{pedigree,TrioSet-method}
\alias{trios,TrioSet-method}
\alias{updateObject,TrioSet-method}
\alias{xyplot,formula,TrioSet-method}

\title{Class \code{"TrioSet"}}
\description{

    A TrioSet is a container for storing high throughput assay data and
  metadata from genotyping arrays when the study design is case-parent
  trios. In our application, de novo copy number alterations in affected
  offspring were of the primary interest.  Examination of the joint
  distribution of the log R ratios and B allele frequencies across
  members in a trio motivates a container with assay data elements that
  are 3-dimensional arrays rather than 2-dimensional matrices.  The
  dimension of the arrays is marker x trio x individual.  Typically, a
  TrioSet instance stores the data for a single chromosome and an
  instance of TrioSetList is a list of TrioSets.  While having a single
  TrioSet for the entire dataset would simplify the classes of objects
  defined in this package, multiple arrays with thousands of trios and
  roughly a million markers are impractical on most machines. In
  addition, storage by chromosome will facilitate parallelization of
  computation that can be carried out independently on different
  chromosomes.

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

\section{Accessors}{
  \describe{
    The \code{object} in the accessor descriptions that follow is a
    \code{TrioSet}:

    \item{}{\code{allNames(object)}: The individual ids for each
    subject. See \code{\link{offspringNames}},
    \code{\link{fatherNames}}, and \code{\link{motherNames}} to list ids
    corresponding to membership in the trio.}

    \item{}{\code{baf(object)}: Extract array of B allele frequencies.}

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

  \item{}{\code{fatherPhenoData(object)}: Extract
    \code{AnnotatedDataFrame} for father pheno data.}



  \item{}{\code{lrr(object) <- value}: assign log R ratios to assayData.
      \code{value} is a 3-dimensional array (feature x trio x sample).}

 \item{}{\code{mindist(object)}: Accessor for the minimum distance
   matrix.}

 \item{}{\code{mindist(object) <- value}: Replacement method for slot
 \code{mindist}.  Value must be a \code{matrix}.  }

 \item{}{\code{motherNames(object) <- value}: Assign a character string
   for the mother identifiers }

  \item{}{\code{motherPhenoData(object)}: Extract
  \code{AnnotatedDataFrame} for mother pheno data.}

   \item{}{\code{motherNames(object)} Return a character vector of
     mother identifiers}

   \item{}{\code{ncol(object)} Number of trios in the \code{TrioSet}
   \code{object}}

   \item{}{\code{offspringNames(object) <- value}: Assign a character
     vector of offspring identifiers}

   \item{}{\code{offspringNames(object)}:  Retrieve character vector of
 offspring identifiers.  Note that the result will identical to
 \code{sampleNames(object)} as the offspring identifiers uniquely
 identify a trio}

  \item{}{\code{offspringPhenoData(object)}: Extract a
  \code{AnnotatedDataFrame} of the sample-level covariates for the
  offspring. Alternatively, use \code{phenoData(object)}. }

   \item{}{\code{order(object)}: order \code{TrioSet} object by
   chromosome and physical position}

  \item{}{\code{pedigree(object)}: accessor for pedigree slot.  Returns
    an object of class \code{Pedigree}.}

  \item{}{\code{trios(object)}: Returns the \code{data.frame} stored in
  slot \code{trios}  of the \code{Pedigree} class object stored in slot
  \code{pedigree} of the \code{TrioSet} object.
}

}
}

\section{Methods for coercing TrioSet objects to another class}{
  \describe{

    \item{}{\code{as(object, "TrioSetList")}: coerce an object of class
    \code{TrioSet} to an object of class \code{TrioSetList}}

    \item{}{\code{as(object, "data.frame")}: coerce an object of class
    \code{TrioSet} to a \code{data.frame} with the following column
    labels: x (position), lrr (log R ratio), gt (genotypes -- ignored if
    genotypes are not included in the \code{assayData}), baf (B allele
    frequencies), id, and is.snp.  The coercion to a \code{data.frame}
    is used primarily for visualization methods using function derived
    from \code{xyplot} in the \pkg{lattice} package.  This is
    experimental and subject to change in future versions.}

  }

}

\section{Compute the minimum distance}{
  \describe{

  \item{}{\code{calculateMindist(object)}: calculate the minimum
    distance for an object of class TrioSet.}

  }

}

\section{Posterior summaries}{

  \describe{

    \item{}{\code{computeBayesFactor(object, ranges, mad.marker,
    mad.sample, returnEmission=FALSE, verbose=TRUE, ...)}: Compute
    posterior probabilities for the trio copy number states. The called
    trio copy number state is the argmax of the posterior
    probabilities. See examples below for computing the standard
    deviation of the markers for argument \code{mad.marker} and the
    standard deviation of the samples for argument \code{mad.sample}.}

  }

}


\section{Visualization}{
  \describe{
    \item{}{\code{xyplot(x, data, ...)}       }
  }
}

\section{Miscellaneous}{
  \describe{

    \item{}{\code{updateObject(object)}: Currently, this method only
    checks the class of the \code{featureData} slot. If the class is
    \code{AnnotatedDataFrame}, the \code{featureData} is updated to the
    \code{GenomeAnnotatedDataFrame} class.}

  }
}

\author{R. Scharpf}


\seealso{
  \code{\link{TrioSet}}, \code{\linkS4class{TrioSetList}}, \code{\link{calculateMindist}}
}

\examples{
	showClass("TrioSet")
	## instantiate a TrioSet with an array of log R ratios (logRR) and B
	## allele frequencies
library(MinimumDistance)
	TrioSet()

        data(trioSetListExample)
        trioSet <- stack(trioSetList)
}
\keyword{classes}
