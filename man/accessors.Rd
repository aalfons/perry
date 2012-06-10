\name{accessors}
\alias{accessors}
\alias{fits}
\alias{fits<-}
\alias{nfits}
\alias{npe}
\alias{peNames}
\alias{peNames<-}
\title{Access or set information on resampling-based prediction error results}
\usage{
  peNames(x)

  peNames(x) <- value

  fits(x)

  fits(x) <- value

  npe(x)

  nfits(x)
}
\arguments{
  \item{x}{an object inheriting from class \code{"perry"}
  or \code{"perrySelect"} that contains prediction error
  results.}

  \item{value}{a vector of replacement values.}
}
\value{
  \code{peNames} returns the names of the prediction error
  results.  The replacement function thereby returns them
  invisibly.

  \code{fits} returns the identifiers of the models for
  objects inheriting from class \code{"perrySelect"} and
  \code{NULL} for objects inheriting from class
  \code{"perry"}.  The replacement function thereby returns
  those values invisibly.

  \code{npe} returns the number of prediction error
  results.

  \code{nfits} returns the number of models included in
  objects inheriting from class \code{"perrySelect"} and
  \code{NULL} for objects inheriting from class
  \code{"perry"}.
}
\description{
  Retrieve or set the names of resampling-based prediction
  error results, retrieve or set the identifiers of the
  models, or retrieve the number of prediction error
  results or included models.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perrySelect}},
  \code{\link{perryTuning}}
}
\keyword{utilities}

