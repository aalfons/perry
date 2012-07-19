\name{subset.perry}
\alias{subset.perry}
\alias{subset.perrySelect}
\title{Subsetting resampling-based prediction error results}
\usage{
  \method{subset}{perry} (x, select = NULL, ...)

  \method{subset}{perrySelect} (x, subset = NULL,
    select = NULL, ...)
}
\arguments{
  \item{x}{an object inheriting from class \code{"perry"}
  or \code{"perrySelect"} that contains prediction error
  results.}

  \item{subset}{a character, integer or logical vector
  indicating the subset of models for which to keep the
  prediction error results.}

  \item{select}{a character, integer or logical vector
  indicating the prediction error results to be extracted.}

  \item{\dots}{currently ignored.}
}
\value{
  An object similar to \code{x} containing just the
  selected results.
}
\description{
  Extract subsets of resampling-based prediction error
  results.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perrySelect}},
  \code{\link{perryTuning}}, \code{\link{subset}}
}
\keyword{utilities}

