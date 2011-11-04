\name{TrioSetList-class}
\Rdversion{1.1}
\docType{class}
\alias{TrioSetList-class}
\alias{allNames,TrioSetList-method}
\alias{allNames}
\alias{[,TrioSetList-method}
\alias{baf,TrioSetList-method}
\alias{baf}
\alias{calculateMindist,TrioSetList-method}
\alias{calculateMindist}
\alias{dims,TrioSetList-method}
\alias{fatherNames,TrioSetList-method}
\alias{fatherNames}
\alias{initialize,TrioSetList-method}
\alias{lrr,TrioSetList-method}
\alias{mad.marker<-,TrioSetList,list-method}
\alias{mad.marker<-,TrioSetList,matrix-method}
%\alias{mad.mindist<-,TrioSetList-method}
\alias{mad.mindist<-,TrioSetList,list-method}
\alias{mad.sample<-,TrioSetList,matrix-method}
\alias{mindist<-,TrioSetList,list-method}
\alias{mindist,TrioSetList-method}
\alias{motherNames,TrioSetList-method}
\alias{motherNames}
\alias{ncol,TrioSetList-method}
\alias{nrow,TrioSetList-method}
\alias{offspringNames,TrioSetList-method}
\alias{offspringNames}
\alias{pedigree,TrioSetList-method}
\alias{pedigree}
\alias{sampleNames,TrioSetList-method}
\alias{sampleSheet,TrioSetList-method}
\alias{sampleSheet}
\alias{stack,TrioSetList-method}
\alias{trios,TrioSetList-method}
\alias{trios}

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

    \item{initialize}{\code{new("TrioSetList")}:  Instantiate an empty
      container.}

    \item{TrioSetList}{ See \code{\link{TrioSetList}} .}

  }

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

\section{Accessors}{

  \describe{

  In the following accessor descriptions, object is a
  \code{TrioSetList}:

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

  \item{}{\code{baf(object)} Returns list of B allele frequencies for each chromosome.}

  \item{}{\code{dims(object)}: Return dimensions of the low-level
  statistical summaries (log R ratios and B allele frequencies) for each
  \code{TrioSet} element in the \code{TrioSetList}.}

  \item{}{\code{fatherNames(object)}: character vector of
    father ids.}

  \item{}{\code{lrr(object)} Returns list of log R ratios for each
    chromosome.}

  \item{}{\code{mad.marker(object) <- value}: The median absolute
  deviation (mad) estimates are assigned to the \code{featureData} of
  each \code{TrioSet} element.  \code{value} can be a \code{list},
  \code{matrix}, or \code{numeric}.  }

  \item{}{\code{mad.mindist(object) <- value}: Standard deviation
  estimates for the minimum distance (across markers for each sample) is
  assigned to each \code{TrioSet} element. \code{value} is a list. Each
  element in the list is the numeric vector of median absolute
  deviations for one chromosome.  }

  \item{}{\code{mad.sample(object) <- value}: Summary measures of the
  standard deviation of the log R ratios for each trio. \code{value} is
  a matrix: rows correspond to trio index, columns correspond to father,
  mother, and offspring.}

  \item{}{\code{mindist(object) <- value}: Assigns the minimum distance
  to each \code{TrioSet} element.  \code{value} is a \code{list}.}

  \item{}{\code{mindist(object)}: Accessor for the minimum distance in
    each \code{TrioSet} element.  Returns a \code{list}.}

  \item{}{\code{motherNames(object)}: Character vector of sample
    identifiers for the mothers.}

  \item{}{\code{ncol(object)}: the number of trios, or equivalently the
    number of offspring.}

  \item{}{\code{nrow(object)}: The number of features across all
  \code{TrioSet} elements.}

   \item{}{\code{offspringNames(object)}: Character vector of offspring
     identifiers. Note that the offspring ids uniquely identify a trio, and the method \code{sampleNames} will
     return the same result.}


   \item{pedigree}{\code{signature(object = "TrioSetList")}: Accessor
     for pedigree information.  See also \code{\linkS4class{Pedigree}}}

   \item{}{\code{sampleNames(object)}: Character vector of unique
   identifiers for father-mother-offspring trio. As the offspring id
   uniquely identifies a trio,\code{sampleNames} returns the same vector
   of ids as \code{offspringNames}}

 \item{}{\code{sampleSheet(object)}: Accessor for metadata of the
   experiment. See \code{\linkS4class{SampleSheet}}.}

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

\section{Computing the minimum distance}{
  \describe{

    The \code{TrioSetList} container is stores the log R ratios for each
    member in the trio.  The 'minimum distance' is the minimum signed
    absolute difference of the parental log R ratios and the offspring
    log R ratios. Specifically, let |O-F| denote the absolute difference
    in the log R ratios comparing offspring to father and |O-M| the
    absolute difference in the log R ratios comparing offspring to
    mother.  The minimum distance at a marker is the signed minimum of
    |O-M| and |O-F|.  Plotted for all markers, deflections from zero can
    indicate a difference in the log R ratio of the offspring and either
    parent (possibly de novo).  For example, a positive minimum distance
    suggests that the log R ratio from the offspring is greater than the
    log R ratio of either parent.  The minimum distance can be computed
    from.
%
\item{}{\code{calculateMindist(object, ..., verbose=TRUE)}:

  The returned value is a \code{list}. Each element in the list contains
  a matrix of the minimum distances from one chromosome (dimension:
  number features x number trios).  Note that the \code{list} can be assigned
  back to the \code{TrioSetList} object using the replacement method
  \code{mindist<-}.  See the example.}

}

}



\section{Combining}{
  \describe{

    A \code{TrioSetList} object contains a list of elements of class
    \code{TrioSet}.  Each list element corresponds to one chromosome. A
    \code{TrioSet} object can be constructed from a \code{TrioSetList}
    by stacking the \code{TrioSet} elements:

    \item{}{\code{stack(object)}: create a \code{TrioSet} object from a
      \code{TrioSetList} object. }

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
data(trioSetListExample)

## accessors
fatherNames(trioSetList)
motherNames(trioSetList)
identical(sampleNames(trioSetList), offspringNames(trioSetList))
pedigree(trioSetList)
sampleSheet(trioSetList)
trios(trioSetList)

## Creating a single TrioSet from the TrioSet elements in the TrioSetList object:
trioSet <- stack(trioSetList)

##
## computing the minimum distance
##
mdlist <- calculateMindist(trioSetList)
mindist(trioSetList) <- mdlist



}
\keyword{classes}
