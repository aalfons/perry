\name{perryTuning}
\alias{perryTuning}
\alias{perryTuning.call}
\alias{perryTuning.function}
\alias{print.perryTuning}
\title{Resampling-based prediction error for tuning parameter selection}
\usage{
  perryTuning(object, ...)

  \method{perryTuning}{function} (object, formula,
    data = NULL, x = NULL, y, tuning = list(),
    args = list(), splits = foldControl(),
    predictFun = predict, predictArgs = list(),
    cost = rmspe, costArgs = list(),
    selectBest = c("min", "hastie"), seFactor = 1,
    names = NULL, envir = parent.frame(), seed = NULL, ...)

  \method{perryTuning}{call} (object, data = NULL,
    x = NULL, y, tuning = list(), splits = foldControl(),
    predictFun = predict, predictArgs = list(),
    cost = rmspe, costArgs = list(),
    selectBest = c("min", "hastie"), seFactor = 1,
    names = NULL, envir = parent.frame(), seed = NULL, ...)
}
\arguments{
  \item{object}{a function or an unevaluated function call
  for fitting a model (see \code{\link{call}} for the
  latter).}

  \item{formula}{a \code{\link[stats]{formula}} describing
  the model.}

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

  \item{tuning}{a list of arguments giving the tuning
  parameter values to be evaluated.  The names of the list
  components should thereby correspond to the argument
  names of the tuning parameters.  For each tuning
  parameter, a vector of values can be supplied.  The
  prediction error is then estimated for all possible
  combinations of tuning parameter values.}

  \item{args}{a list of additional arguments to be passed
  to the model fitting function.}

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

  \item{predictFun}{a function to compute predictions for
  the test data.  It should expect the fitted model to be
  passed as the first argument and the test data as the
  second argument, and must return either a vector or a
  matrix containing the predicted values.  The default is
  to use the \code{\link[stats]{predict}} method of the
  fitted model.}

  \item{predictArgs}{a list of additional arguments to be
  passed to \code{predictFun}.}

  \item{cost}{a cost function measuring prediction loss.
  It should expect the observed values of the response to
  be passed as the first argument and the predicted values
  as the second argument, and must return either a
  non-negative scalar value, or a list with the first
  component containing the prediction error and the second
  component containing the standard error.  The default is
  to use the root mean squared prediction error (see
  \code{\link{cost}}).}

  \item{costArgs}{a list of additional arguments to be
  passed to the prediction loss function \code{cost}.}

  \item{selectBest}{a character string specifying a
  criterion for selecting the best model.  Possible values
  are \code{"min"} (the default) or \code{"hastie"}.  The
  former selects the model with the smallest prediction
  error.  The latter is useful for models with a tuning
  parameter controlling the complexity of the model (e.g.,
  penalized regression).  It selects the most parsimonious
  model whose prediction error is no larger than
  \code{seFactor} standard errors above the prediction
  error of the best overall model.  Note that the models
  are thereby assumed to be ordered from the most
  parsimonious one to the most complex one.  In particular
  a one-standard-error rule is frequently applied.}

  \item{seFactor}{a numeric value giving a multiplication
  factor of the standard error for the selection of the
  best model.  This is ignored if \code{selectBest} is
  \code{"min"}.}

  \item{names}{an optional character vector giving names
  for the arguments containing the data to be used in the
  function call (see \dQuote{Details}).}

  \item{envir}{the \code{\link{environment}} in which to
  evaluate the function call for fitting the models (see
  \code{\link{eval}}).}

  \item{seed}{optional initial seed for the random number
  generator (see \code{\link{.Random.seed}}).}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  If \code{tuning} is an empty list, \code{\link{perryFit}}
  is called to return an object of class \code{"perry"}.

  Otherwise an object of class \code{"perryTuning"} (which
  inherits from class \code{"perrySelect"}) with the
  following components is returned:

  \item{pe}{a data frame containing the estimated
  prediction errors for all combinations of tuning
  parameter values.  In case of more than one replication,
  those are average values over all replications.}

  \item{se}{a data frame containing the estimated standard
  errors of the prediction loss for all combinations of
  tuning parameter values.}

  \item{reps}{a data frame containing the estimated
  prediction errors from all replications for all
  combinations of tuning parameter values.  This is only
  returned in case of more than one replication.}

  \item{splits}{an object giving the data splits used to
  estimate the prediction error.}

  \item{y}{the response.}

  \item{yHat}{a list containing the predicted values for
  all combinations of tuning parameter values.  Each list
  component is again a list containing the corresponding
  predicted values from all replications.}

  \item{best}{an integer vector giving the indices of the
  optimal combinations of tuning parameters.}

  \item{selectBest}{a character string specifying the
  criterion used for selecting the best model.}

  \item{seFactor}{a numeric value giving the multiplication
  factor of the standard error used for the selection of
  the best model.}

  \item{tuning}{a data frame containing the grid of tuning
  parameter values for which the prediction error was
  estimated.}

  \item{seed}{the seed of the random number generator
  before estimation of the prediction error.}

  \item{call}{the matched function call.}
}
\description{
  Select tuning parameters of a model by estimating the
  respective prediction errors via (repeated) \eqn{K}-fold
  cross-validation, (repeated) random splitting (also known
  as random subsampling or Monte Carlo cross-validation),
  or the bootstrap.  It is thereby possible to supply a
  model fitting function or an unevaluated function call to
  a model fitting function.
}
\details{
  (Repeated) \eqn{K}-fold cross-validation is performed in
  the following way.  The data are first split into \eqn{K}
  previously obtained blocks of approximately equal size
  (given by \code{folds}).  Each of the \eqn{K} data blocks
  is left out once to fit the model, and predictions are
  computed for the observations in the left-out block with
  \code{predictFun}.  Thus a prediction is obtained for
  each observation.  The response variable and the obtained
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

  In any case, if the response is a vector but
  \code{predictFun} returns a matrix, the prediction error
  is computed for each column.  A typical use case for this
  behavior would be if \code{predictFun} returns
  predictions from an initial model fit and stepwise
  improvements thereof.

  If \code{formula} or \code{data} are supplied, all
  variables required for fitting the models are added as
  one argument to the function call, which is the typical
  behavior of model fitting functions with a
  \code{\link[stats]{formula}} interface.  In this case,
  the accepted values for \code{names} depend on the
  method.  For the \code{function} method, a character
  vector of length two should supplied, with the first
  element specifying the argument name for the formula and
  the second element specifying the argument name for the
  data (the default is to use \code{c("formula", "data")}).
  Note that names for both arguments should be supplied
  even if only one is actually used.  For the \code{call}
  method, which does not have a \code{formula} argument, a
  character string specifying the argument name for the
  data should be supplied (the default is to use
  \code{"data"}).

  If \code{x} is supplied, on the other hand, the predictor
  matrix and the response are added as separate arguments
  to the function call.  In this case, \code{names} should
  be a character vector of length two, with the first
  element specifying the argument name for the predictor
  matrix and the second element specifying the argument
  name for the response (the default is to use \code{c("x",
  "y")}).  It should be noted that the \code{formula} or
  \code{data} arguments take precedence over \code{x}.
}
\note{
  The same data splits are used for all combinations of
  tuning parameter values for maximum comparability.
}
\author{
  Andreas Alfons
}
\references{
  Hastie, T., Tibshirani, R. and Friedman, J. (2009)
  \emph{The Elements of Statistical Learning: Data Mining,
  Inference, and Prediction}.  Springer, 2nd edition.
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perrySelect}},
  \code{\link{cvFolds}}, \code{\link{randomSplits}},
  \code{\link{bootSamples}}, \code{\link{cost}}
}
\keyword{utilities}

