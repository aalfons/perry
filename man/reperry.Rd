\name{reperry}
\alias{reperry}
\alias{reperry.perry}
\alias{reperry.perrySelect}
\title{Recompute resampling-based prediction error measures}
\usage{
  reperry(object, ...)

  \method{reperry}{perry} (object, cost = rmspe, ...)

  \method{reperry}{perrySelect} (object, cost = rmspe, ...)
}
\arguments{
  \item{object}{an object inheriting from class
  \code{"perry"} or \code{"perrySelect"} that contains
  prediction error results.}

  \item{cost}{a cost function measuring prediction loss.
  It should expect the observed values of the response to
  be passed as the first argument and the predicted values
  as the second argument, and must return either a
  non-negative scalar value, or a list with the first
  component containing the prediction error and the second
  component containing the standard error.  The default is
  to use the root mean squared prediction error (see
  \code{\link{cost}}).}

  \item{\dots}{for the generic function, additional
  arguments to be passed down to methods.  For the methods,
  additional arguments to be passed to the prediction loss
  function \code{cost}.}
}
\value{
  An object similar to \code{object} containing the results
  for the new measure of prediction loss.
}
\description{
  Recompute prediction error measures for previously
  obtained objects that contain resampling-based prediction
  error results.  This is useful for computing a different
  measure of prediction loss.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perryTuning}},
  \code{\link{perrySelect}}
}
\keyword{utilities}

