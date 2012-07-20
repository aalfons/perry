\name{perry}
\alias{perry}
\alias{perry.lm}
\alias{perry.lmrob}
\alias{perry.lts}
\title{Resampling-based prediction error for fitted models}
\usage{
  perry(object, ...)

  \method{perry}{lm} (object, splits = foldControl(),
    cost = rmspe, seed = NULL, ...)

  \method{perry}{lmrob} (object, splits = foldControl(),
    cost = rtmspe, seed = NULL, ...)

  \method{perry}{lts} (object, splits = foldControl(),
    fit = c("reweighted", "raw", "both"), cost = rtmspe,
    seed = NULL, ...)
}
\arguments{
  \item{object}{the fitted model for which to estimate the
  prediction error.}

  \item{splits}{an object of class \code{"cvFolds"} (as
  returned by \code{\link{cvFolds}}) or a control object of
  class \code{"foldControl"} (see
  \code{\link{foldControl}}) defining the folds of the data
  for (repeated) \eqn{K}-fold cross-validation, an object
  of class \code{"randomSplits"} (as returned by
  \code{\link{randomSplits}}) or a control object of class
  \code{"splitControl"} (see \code{\link{splitControl}})
  defining random data splits, or an object of class
  \code{"bootSamples"} (as returned by
  \code{\link{bootSamples}}) or a control object of class
  \code{"bootControl"} (see \code{\link{bootControl}})
  defining bootstrap samples.}

  \item{fit}{a character string specifying for which fit to
  estimate the prediction error.  Possible values are
  \code{"reweighted"} (the default) for the prediction
  error of the reweighted fit, \code{"raw"} for the
  prediction error of the raw fit, or \code{"both"} for the
  prediction error of both fits.}

  \item{cost}{a cost function measuring prediction loss.
  It should expect the observed values of the response to
  be passed as the first argument and the predicted values
  as the second argument, and must return either a
  non-negative scalar value, or a list with the first
  component containing the prediction error and the second
  component containing the standard error.  The default is
  to use the root mean squared prediction error for the
  \code{"lm"} method and the root trimmed mean squared
  prediction error for the \code{"lmrob"} and \code{"lts"}
  methods (see \code{\link{cost}}).}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}

  \item{\dots}{for the generic function, additional
  arguments to be passed down to methods.  For the methods,
  additional arguments to be passed to the prediction loss
  function \code{cost}.}
}
\value{
  An object of class \code{"perry"} with the following
  components:

  \item{pe}{a numeric vector containing the estimated
  prediction errors.  For the \code{"lm"} and
  \code{"lmrob"} methods, this is a single numeric value.
  For the \code{"lts"} method, this contains one value for
  each of the requested fits.  In case of more than one
  replication, those are average values over all
  replications.}

  \item{se}{a numeric vector containing the estimated
  standard errors of the prediction loss.  For the
  \code{"lm"} and \code{"lmrob"} methods, this is a single
  numeric value.  For the \code{"lts"} method, this
  contains one value for each of the requested fits.}

  \item{reps}{a numeric matrix containing the estimated
  prediction errors from all replications.  For the
  \code{"lm"} and \code{"lmrob"} methods, this is a matrix
  with one column.  For the \code{"lts"} method, this
  contains one column for each of the requested fits.
  However, this is only returned in case of more than one
  replication.}

  \item{splits}{an object giving the data splits used to
  estimate the prediction error.}

  \item{y}{the response.}

  \item{yHat}{a list containing the predicted values from
  all replications.}

  \item{seed}{the seed of the random number generator
  before estimation of the prediction error.}

  \item{call}{the matched function call.}
}
\description{
  Estimate the prediction error of a fitted model via
  (repeated) \eqn{K}-fold cross-validation, (repeated)
  random splitting (also known as random subsampling or
  Monte Carlo cross-validation), or the bootstrap.  Methods
  are available for least squares fits computed with
  \code{\link[stats]{lm}} as well as for the following
  robust alternatives: MM-type models computed with
  \code{\link[robustbase]{lmrob}} and least trimmed squares
  fits computed with \code{\link[robustbase]{ltsReg}}.
}
\note{
  The \code{perry} methods extract the data from the fitted
  model and call \code{\link{perryFit}} to perform
  resampling-based prediction error estimation.  Users may
  prefer the wrapper functions \code{\link{repCV}},
  \code{\link{repRS}} and \code{\link{bootPE}}.
}
\examples{
## load data and fit an LS regression model
data("mtcars")
fit <- lm(mpg ~ wt + cyl, data=mtcars)

## perform cross-validation
# K-fold CV
perry(fit, foldControl(K = 5, R = 10), seed = 1234)
# leave-one-out CV
perry(fit, foldControl(K = nrow(mtcars)))

## perform random splitting
perry(fit, splitControl(m = 6, R = 10), seed = 1234)

## perform bootstrap prediction error estimation
# 0.632 estimator
perry(fit, bootControl(R = 10, type = "0.632"), seed = 1234)
# out-of-bag estimator
perry(fit, bootControl(R = 10, type = "out-of-bag"), seed = 1234)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{repCV}},
  \code{\link{repRS}}, \code{\link{bootPE}}
}
\keyword{utilities}
