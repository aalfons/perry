\name{perryTool}
\alias{perryTool}
\title{Low-level function for resampling-based prediction error measures}
\usage{
  perryTool(call, data = NULL, x = NULL, y, cost = rmspe,
    splits, names = NULL, predictArgs = list(),
    costArgs = list(), envir = parent.frame())
}
\arguments{
  \item{call}{an unevaluated function call for fitting a
  model (see \code{\link{call}}).}

  \item{data}{a data frame containing the variables
  required for fitting the models.  This is typically used
  if the model in the function call is described by a
  \code{\link[stats]{formula}}.}

  \item{x}{a numeric matrix containing the predictor
  variables.  This is typically used if the function call
  for fitting the models requires the predictor matrix and
  the response to be supplied as separate arguments.}

  \item{y}{a numeric vector or matrix containing the
  response.}

  \item{cost}{a cost function measuring prediction loss.
  It should expect the observed values of the response to
  be passed as the first argument and the predicted values
  as the second argument, and must return either a
  non-negative scalar value, or a list with the first
  component containing the prediction error and the second
  component containing the standard error.  The default is
  to use the root mean squared prediction error (see
  \code{\link{cost}}).}

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

  \item{names}{an optional character vector giving names
  for the arguments containing the data to be used in the
  function call (see \dQuote{Details}).}

  \item{predictArgs}{a list of additional arguments to be
  passed to the \code{\link[stats]{predict}} method of the
  fitted models.}

  \item{costArgs}{a list of additional arguments to be
  passed to the prediction loss function \code{cost}.}

  \item{envir}{the \code{\link{environment}} in which to
  evaluate the function call for fitting the models (see
  \code{\link{eval}}).}
}
\value{
  If only one replication is requested and the prediction
  loss function \code{cost} also returns the standard
  error, a list is returned, with the first component
  containing the estimated prediction errors and the second
  component the corresponding estimated standard errors.

  Otherwise the return value is a numeric matrix in which
  each column contains the respective estimated prediction
  errors from all replications.
}
\description{
  Basic function to estimate the prediction error of a
  model via (repeated) \eqn{K}-fold cross-validation,
  (repeated) random splitting (also known as random
  subsampling or Monte Carlo cross-validation), or the
  bootstrap.  The model is thereby specified by an
  unevaluated function call to a model fitting function.
}
\details{
  (Repeated) \eqn{K}-fold cross-validation is performed in
  the following way.  The data are first split into \eqn{K}
  previously obtained blocks of approximately equal size
  (given by \code{folds}).  Each of the \eqn{K} data blocks
  is left out once to fit the model, and predictions are
  computed for the observations in the left-out block with
  the \code{\link[stats]{predict}} method of the fitted
  model.  Thus a prediction is obtained for each
  observation.  The response variable and the obtained
  predictions for all observations are then passed to the
  prediction loss function \code{cost} to estimate the
  prediction error.  For repeated \eqn{K}-fold
  cross-validation (as indicated by \code{splits}), this
  process is replicated and the estimated prediction errors
  from all replications are returned.

  (Repeated) random splitting is performed similarly.  In
  each replication, the data are split into a training set
  and a test set at random.  Then the training data is used
  to fit the model, and predictions are computed for the
  test data.  Hence only the response values from the test
  data and the corresponding predictions are passed to the
  prediction loss function \code{cost}.

  For the bootstrap estimator, each bootstrap sample is
  used as training data to fit the model.  The out-of-bag
  estimator uses the observations that do not enter the
  bootstrap sample as test data and computes the prediction
  loss function \code{cost} for those out-of-bag
  observations.  The 0.632 estimator is computed as a
  linear combination of the out-of-bag estimator and the
  prediction loss of the fitted values of the model
  computed from the full sample.

  In any case, if the response is a vector but the
  \code{\link[stats]{predict}} method of the fitted models
  returns a matrix, the prediction error is computed for
  each column.  A typical use case for this behavior would
  be if the \code{\link[stats]{predict}} method returns
  predictions from an initial model fit and stepwise
  improvements thereof.

  If \code{data} is supplied, all variables required for
  fitting the models are added as one argument to the
  function call, which is the typical behavior of model
  fitting functions with a \code{\link[stats]{formula}}
  interface.  In this case, a character string specifying
  the argument name can be passed via \code{names} (the
  default is to use \code{"data"}).

  If \code{x} is supplied, on the other hand, the predictor
  matrix and the response are added as separate arguments
  to the function call.  In this case, \code{names} should
  be a character vector of length two, with the first
  element specifying the argument name for the predictor
  matrix and the second element specifying the argument
  name for the response (the default is to use \code{c("x",
  "y")}).  It should be noted that \code{data} takes
  precedence over \code{x} if both are supplied.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perryTuning}},
  \code{\link{cvFolds}}, \code{\link{foldControl}},
  \code{\link{randomSplits}}, \code{\link{splitControl}},
  \code{\link{bootSamples}}, \code{\link{bootControl}},
  \code{\link{cost}}
}
\keyword{utilities}

