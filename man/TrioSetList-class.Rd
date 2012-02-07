\name{TrioSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{TrioSetList-class}
\alias{allNames,TrioSetList-method}
\alias{allNames}
\alias{[,TrioSetList-method}
\alias{$,TrioSetList-method}
\alias{[[,TrioSetList-method}
\alias{annotation,TrioSetList-method}
\alias{baf,TrioSetList-method}
\alias{baf}
\alias{computeBayesFactor,TrioSetList-method}
\alias{dims,TrioSetList-method}
\alias{fatherNames,TrioSetList-method}
\alias{fatherPhenoData,TrioSetList-method}
\alias{motherPhenoData,TrioSetList-method}
\alias{offspringPhenoData,TrioSetList-method}
\alias{fatherNames}
\alias{initialize,TrioSetList-method}
\alias{length,TrioSetList-method}
\alias{lrr,TrioSetList-method}
\alias{mindist<-,TrioSetList,list-method}
\alias{mindist,TrioSetList-method}
\alias{motherNames,TrioSetList-method}
\alias{motherNames}
\alias{ncol,TrioSetList-method}
\alias{nrow,TrioSetList-method}
\alias{offspringNames,TrioSetList-method}
\alias{offspringNames}
\alias{order,TrioSetList-method}
\alias{pedigree,TrioSetList-method}
\alias{pedigree}
\alias{phenoData,TrioSetList-method}
\alias{sampleNames,TrioSetList-method}
\alias{stack,TrioSetList-method}
\alias{trios,TrioSetList-method}
\alias{trios}
\alias{updateObject,TrioSetList-method}

\title{Class \code{"TrioSetList"}}
\description{
  A container storing pedigree information, as well as low-level
  statistical summaries used for copy number estimation: the log R
  ratios and B allele frequencies.  The list structure is organized by
  chromosome, where each element of the list is a \code{TrioSet} object.
}


\section{Objects from the Class}{

  \describe{

    Objects from the class can be initialized by:

    \item{}{\code{new("TrioSetList")}:  Instantiate an empty
      container.}

    \item{}{\code{TrioSetList()}: See \code{\link{TrioSetList}} .}
  }
}

\section{Slots}{
  \describe{
    \item{\code{assayDataList}:}{list of arrays containing log R ratios and B
      allele frequencies. }
    \item{\code{featureDataList}:}{list of containing feature
  annotation. Each element of the list is a \code{AnnotatedDataFrame}.}
    \item{\code{pedigree}:}{ Object of class \code{"Pedigree"}.  Contains
      information on the trio-relationships.}
    \item{\code{motherPhenoData}:}{ Object of class \code{"AnnotatedDataFrame"}.  Contains
      sample-level covariates for the mother.}
    \item{\code{fatherPhenoData}:}{ Object of class \code{"AnnotatedDataFrame"}.  Contains
      sample-level covariates for the father.}
    \item{\code{phenoData}:}{ Object of class \code{"AnnotatedDataFrame"}.  Contains
      sample-level covariates for the offspring.}
    \item{\code{chromosome}:}{ Integer vector indicating which autosomes
      are contained in the \code{TrioSetList} object.}
  }
}

\section{Accessors}{

  \describe{

  In the following accessor descriptions, object is a
  \code{TrioSetList}:

  \item{}{\code{object$NAME}: Extract phenotype 'NAME' for offspring.}

  \item{}{\code{"[["}: Extract a \code{TrioSet} object.}

  \item{}{\code{allnames(object)}: Returns character vector of
 all the sample names.  Note that \code{sampleNames} of a
 \code{TrioSetList} object is not the same.  In particular,
 \code{sampleNames(object)} returns only the character vector of
 offspring ids which uniquely identify a trio.  Hence, a separate
 method, \code{allNames}, is supplied when all the sample ids in the
 data set are required.  Finally, note that the length of the vector
 returned by \code{allNames} is not necessarily a multiple of 3 as
 mothers and fathers with multiple offspring would be included in
 multiple trios.}

\item{}{\code{annotation(object)}: character string indicating the
  array platform }

  \item{}{\code{baf(object)}: Returns a list of B allele frequencies for
    each chromosome. Each element in the list is a 3-dimensional array
    (features x trios x samples).}

  \item{}{\code{dims(object)}: Return dimensions of the low-level
  statistical summaries (log R ratios and B allele frequencies) for each
  \code{TrioSet} element in the \code{TrioSetList}.}

  \item{}{\code{fatherNames(object)}: character vector of
    father ids.}

  \item{}{\code{fatherPhenoData(object)}: Extract a
  \code{AnnotatedDataFrame} of the sample-level covariates for the
  father.}

  \item{}{\code{length(object)}: The number of chromosomes.}

  \item{}{\code{lrr(object)} Returns list of log R ratios for each
    chromosome.  Each element in the list is a 3-dimensional array
    (features x trios x samples).}

  \item{}{\code{mindist(object) <- value}: Assigns the minimum distance
  to each \code{TrioSet} element.  \code{value} is a \code{list}.}

  \item{}{\code{mindist(object)}: Accessor for the minimum distance in
    each \code{TrioSet} element.  Returns a \code{list}.}

  \item{}{\code{motherNames(object)}: Character vector of sample
    identifiers for the mothers.}

  \item{}{\code{motherPhenoData(object)}: Extract a
  \code{AnnotatedDataFrame} of the sample-level covariates for the
  mother.}

  \item{}{\code{ncol(object)}: the number of trios, or equivalently the
    number of offspring.}

  \item{}{\code{nrow(object)}: The number of features across all
  \code{TrioSet} elements.}

   \item{}{\code{offspringNames(object)}: Character vector of offspring
     identifiers. Note that the offspring ids uniquely identify a trio, and the method \code{sampleNames} will
     return the same result.}

   \item{}{\code{order(object)}: order \code{TrioSet} object by
   chromosome and physical position}


   \item{pedigree}{\code{signature(object = "TrioSetList")}: Accessor
     for pedigree information.  See also \code{\linkS4class{Pedigree}}}

  \item{}{\code{offspringPhenoData(object)}: Extract a
  \code{AnnotatedDataFrame} of the sample-level covariates for the
  offspring. Alternatively, use \code{phenoData(object)}. }

  \item{}{\code{phenoData(object)}: Extract a
  \code{AnnotatedDataFrame} of the sample-level covariates for the offspring.}

   \item{}{\code{sampleNames(object)}: Character vector of unique
   identifiers for father-mother-offspring trio. As the offspring id
   uniquely identifies a trio,\code{sampleNames} returns the same vector
   of ids as \code{offspringNames}}

% \item{}{\code{sampleSheet(object)}: Accessor for metadata of the
%   experiment. See \code{\linkS4class{SampleSheet}}.}

 \item{}{\code{trios(object)}: Returns a
 \code{data.frame} of the trios. Each row in the {data.frame} contains
 the sample identifiers for the father, mother, and offspring.  Parents
 with multiple offspring will appear in multiple rows.}

}

}

\section{Subsetting}{

  \describe{

    \item{}{ \code{x[i, j]}: \code{i} selects the list elements.
    \code{j} selects the trio for each list element.  The list elements
    have class \code{\linkS4class{TrioSet}}.}

  }

}



%\section{Computing the minimum distance}{
%  \describe{
%
%    The \code{TrioSetList} container stores the log R ratios and B
%    allele frequencies for each member in the trio.
%%
%\item{}{\code{calculateMindist(object, ..., verbose=TRUE)}:
%
%  The returned value is a \code{list}. Each element in the list contains
%  a matrix of the minimum distances from one chromosome (dimension:
%  number features x number trios).  Note that the \code{list} can be assigned
%  back to the \code{TrioSetList} object using the replacement method
%  \code{mindist<-}.  See the example.}
%
%}
%
%}



\section{Stacking a \code{TrioSetList} object}{
  \describe{

    A \code{TrioSetList} object contains a list of elements of class
    \code{TrioSet}.  Each list element corresponds to one chromosome. A
    \code{TrioSet} object can be constructed from a \code{TrioSetList}
    by stacking the \code{TrioSet} elements:

    \item{}{\code{stack(object)}: creates a \code{TrioSet} object from a
      \code{TrioSetList} object. }
  }
}

\section{Compute the minimum distance}{
  \describe{
    \item{}{\code{calculateMindist(object)}: Compute the minimum distance
  from the list of log R ratio arrays.  Returns a list of matrices of
  the minimum distance; each element in the list is the minimum distance
  for one chromosome. }

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

\section{Miscellaneous}{
  \describe{

    \item{}{\code{updateObject(object)}: Currently, this method only
    checks the class of the elements in the \code{featureDataList}
    slot.  If the elements are \code{AnnotatedDataFrame} instead of
    \code{GenomeAnnotatedDataFrame}, the \code{featureDataList} slot is updated.}

  }
}


\author{
R. Scharpf
}

\seealso{
  \code{\linkS4class{TrioSet}}, \code{\linkS4class{Pedigree}}
}

\examples{
showClass("TrioSetList")
TrioSetList()
data(trioSetExample)
}
\keyword{classes}
