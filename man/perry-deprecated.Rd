\name{perry-deprecated}
\alias{bootPE}
\alias{perry-deprecated}
\alias{repCV}
\alias{repRS}
\title{Deprecated functions in package \pkg{perry}}
\usage{
  repCV(object, K = 5, R = 1,
    foldType = c("random", "consecutive", "interleaved"),
    grouping = NULL, folds = NULL, ...)

  repRS(object, m, R = 1, grouping = NULL, splits = NULL,
    ...)

  bootPE(object, R = 1,
    bootType = c("0.632", "out-of-bag"), grouping = NULL,
    samples = NULL, ...)
}
\arguments{
  \item{object}{the fitted model for which to estimate the
  prediction error.}

  \item{K}{an integer giving the number of folds into which
  the observations should be split (the default is five).
  Setting \code{K} equal to the number of observations or
  groups yields leave-one-out cross-validation.}

  \item{m}{an integer giving the number of observations or
  groups of observations to be used as test data.}

  \item{R}{an integer giving the number of replications.
  In \code{repCV}, this is ignored for for leave-one-out
  cross-validation and other non-random splits of the
  data.}

  \item{foldType}{a character string specifying the type of
  folds to be generated.  Possible values are
  \code{"random"} (the default), \code{"consecutive"} or
  \code{"interleaved"}.}

  \item{bootType}{a character string specifying a bootstrap
  estimator.  Possible values are \code{"0.632"} (the
  default), or \code{"out-of-bag"}.}

  \item{grouping}{a factor specifying groups of
  observations.  If supplied, the data are split according
  to the groups rather than individual observations such
  that all observations within a group belong either to the
  training or test data.}

  \item{folds}{an object of class \code{"cvFolds"} (as
  returned by \code{\link{cvFolds}}) or a control object of
  class \code{"foldControl"} (see
  \code{\link{foldControl}}) defining the folds of the data
  for (repeated) \eqn{K}-fold cross-validation.  If
  supplied, this is preferred over the arguments for
  generating cross-validation folds.}

  \item{splits}{an object of class \code{"randomSplits"}
  (as returned by \code{\link{randomSplits}}) or a control
  object of class \code{"splitControl"} (see
  \code{\link{splitControl}}) defining random data splits.
  If supplied, this is preferred over the arguments for
  generating random data splits.}

  \item{samples}{an object of class \code{"bootSamples"}
  (as returned by \code{\link{bootSamples}}) or a control
  object of class \code{"bootControl"} (see
  \code{\link{bootControl}}) defining bootstrap samples.
  If supplied, this is preferred over the arguments for
  generating bootstrap samples.}

  \item{\dots}{additional arguments to be passed down to
  \code{\link{perry}}.}
}
\description{
  These functions are provided for compatibility with older
  versions only, and may be defunct as soon as the next
  release.
}
\details{
  \code{repCV}, \code{repRS} and \code{bootPE} are wrapper
  functions for \code{\link{perry}} that perform (repeated)
  \eqn{K}-fold cross-validation, (repeated) random
  splitting (also known as random subsampling or Monte
  Carlo cross-validation) and the bootstrap, respectively.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[base]{Deprecated}}
}
\keyword{utilities}

