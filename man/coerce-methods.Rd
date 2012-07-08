\name{coerceMethods}
\alias{coerce,RangedDataCNV,GRanges-method}

\title{Methods for coercing objects from one class to another}

\description{

  Methods for coercing objects between classes.

}

\section{Methods for RangedDataCNV}{

  \code{RangedDataCNV}-derived classes will eventually be deprecated in
  \pkg{oligoClasses} and replaced by \code{GRanges}.  To coerce from a
  \code{RangedDataCNV} object to a \code{GRangs} object, use

  \describe{

    \item{}{
      \code{as(object, "GRanges")}:

      The coercion to GRanges retains information on coverage
      ('numberProbes'), sample identifiers ('sample'), state (when
      applicable), and segment means (when applicable).

    }
  }

}

\author{R. Scharpf}

\keyword{methods}