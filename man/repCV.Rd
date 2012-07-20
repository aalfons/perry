\name{repCV}
\alias{repCV}
\title{(Repeated) cross-validation for fitted models}
\usage{
  repCV(object, K = 5, R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    grouping = NULL, folds = NULL, ...)
}
\arguments{
  \item{object}{the fitted model for which to estimate the
  prediction error.}

  \item{K}{an integer giving the number of folds into which
  the observations should be split (the default is five).
  Setting \code{K} equal to the number of observations or
  groups yields leave-one-out cross-validation.}

  \item{R}{an integer giving the number of replications for
  repeated \eqn{K}-fold cross-validation.  This is ignored
  for for leave-one-out cross-validation and other
  non-random splits of the data.}

  \item{foldType}{a character string specifying the type of
  folds to be generated.  Possible values are
  \code{"random"} (the default), \code{"consecutive"} or
  \code{"interleaved"}.}

  \item{grouping}{a factor specifying groups of
  observations.  If supplied, the data are split according
  to the groups rather than individual observations such
  that all observations within a group belong to the same
  fold.}

  \item{folds}{an object of class \code{"cvFolds"} (as
  returned by \code{\link{cvFolds}}) or a control object of
  class \code{"foldControl"} (see
  \code{\link{foldControl}}) defining the folds of the data
  for (repeated) \eqn{K}-fold cross-validation.  If
  supplied, this is preferred over the arguments for
  generating cross-validation folds.}

  \item{\dots}{additional arguments to be passed down to
  \code{\link{perry}}.}
}
\value{
  An object of class \code{"perry"} as returned by
  \code{\link{perry}}.
}
\description{
  Estimate the prediction error of a fitted model via
  (repeated) \eqn{K}-fold cross-validation.  This works for
  any model for which a \code{\link{perry}} method is
  available.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perry}}, \code{\link{repRS}},
  \code{\link{bootPE}}
}
\keyword{utilities}

