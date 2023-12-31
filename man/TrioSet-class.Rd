% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllClasses.R, R/methods-TrioSet.R,
%   R/methods-ff_array.R
\docType{class}
\name{TrioSet-class}
\alias{TrioSet-class}
\alias{pedigree,TrioSet-method}
\alias{show,TrioSet-method}
\alias{mindist,TrioSet-method}
\alias{mindist<-,TrioSet,matrix-method}
\alias{dim,TrioSet-method}
\alias{trios,TrioSet-method}
\alias{[,TrioSet,ANY,ANY,ANY-method}
\alias{"[",TrioSet,ANY-method}
\alias{MAP,TrioSet,GRanges-method}
\alias{[,ff_array,ANY,ANY,ANY-method}
\alias{"[",ff_array,ANY-method}
\title{Deprecated class for storing low-level genomic data for trios}
\usage{
\S4method{pedigree}{TrioSet}(object)

\S4method{show}{TrioSet}(object)

\S4method{mindist}{TrioSet}(object)

\S4method{mindist}{TrioSet,matrix}(object) <- value

\S4method{dim}{TrioSet}(x)

\S4method{trios}{TrioSet}(object)

\S4method{[}{TrioSet,ANY,ANY,ANY}(x, i, j, ..., drop = FALSE)

\S4method{MAP}{TrioSet,GRanges}(
  object,
  ranges,
  transition_param = TransitionParam(),
  emission_param = EmissionParam(),
  mdThr = 0.9,
  ...
)

\S4method{[}{ff_array,ANY,ANY,ANY}(x, i, j, ..., drop = FALSE)
}
\arguments{
\item{object}{a \code{TrioSet} object}

\item{value}{a \code{matrix}}

\item{x}{a \code{TrioSet} object}

\item{i}{a numeric vector for subsetting rows  (optional)}

\item{j}{a numeric vector for subsetting trios (optional)}

\item{...}{additional arguments passed to subsetting methods for matrices and data frames}

\item{drop}{logical. Whether to simplify matrices to numeric
vectors.  This should be left as FALSE.}

\item{ranges}{a \code{GRanges} object}

\item{transition_param}{an object of class \code{TransitionParam}}

\item{emission_param}{an object of class \code{EmissionParam}}

\item{mdThr}{the minimum absolute value of the minimum distance
segment mean. Segments with means below \code{mdThr} in absolute
value will not be called as they are unlikely to be de novo.}
}
\description{
This class is deprecated and will be defunct in a future release.
}
\section{Slots}{

\describe{
\item{\code{fatherPhenoData}}{\code{AnnotatedDataFrame} containing covariates for the father}

\item{\code{motherPhenoData}}{\code{AnnotatedDataFrame} containing covariates for the mother}

\item{\code{pedigree}}{an object of class \code{Pedigree}}

\item{\code{mindist}}{a numeric matrix of the minimum distance for each trio, or NULL}
}}

