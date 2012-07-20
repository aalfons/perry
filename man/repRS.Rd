\name{repRS}
\alias{repRS}
\title{(Repeated) random splitting for fitted models}
\usage{
  repRS(object, m, R = 1, grouping = NULL, splits = NULL,
    ...)
}
\arguments{
  \item{object}{the fitted model for which to estimate the
  prediction error.}

  \item{m}{an integer giving the number of observations or
  groups of observations to be used as test data.}

  \item{R}{an integer giving the number of random data
  splits.}

  \item{grouping}{a factor specifying groups of
  observations.  If supplied, the data are split according
  to the groups rather than individual observations such
  that all observations within a group belong either to the
  training or test data.}

  \item{splits}{an object of class \code{"randomSplits"}
  (as returned by \code{\link{randomSplits}}) or a control
  object of class \code{"splitControl"} (see
  \code{\link{splitControl}}) defining random data splits.
  If supplied, this is preferred over the arguments for
  generating random data splits.}

  \item{\dots}{additional arguments to be passed down to
  \code{\link{perry}}.}
}
\value{
  An object of class \code{"perry"} as returned by
  \code{\link{perry}}.
}
\description{
  Estimate the prediction error of a fitted model via
  (repeated) random splitting (also known as random
  subsampling or Monte Carlo cross-validation).  This works
  for any model for which a \code{\link{perry}} method is
  available.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perry}}, \code{\link{repCV}},
  \code{\link{bootPE}}
}
\keyword{utilities}

