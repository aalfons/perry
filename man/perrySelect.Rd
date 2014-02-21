\name{perrySelect}
\alias{perrySelect}
\alias{print.perrySelect}
\title{Model selection via resampling-based prediction error measures}
\usage{
  perrySelect(..., .list = list(...), .reshape = FALSE,
    .selectBest = c("min", "hastie"), .seFactor = 1)
}
\arguments{
  \item{\dots}{objects inheriting from class \code{"perry"}
  or \code{"perrySelect"} that contain prediction error
  results.}

  \item{.list}{a list of objects inheriting from class
  \code{"perry"} or \code{"perrySelect"}.  If supplied,
  this is preferred over objects supplied via the \dots
  argument.}

  \item{.reshape}{a logical indicating whether objects with
  more than one column of prediction error results should
  be reshaped to have only one column (see
  \dQuote{Details}).}

  \item{.selectBest}{a character string specifying a
  criterion for selecting the best model.  Possible values
  are \code{"min"} (the default) or \code{"hastie"}.  The
  former selects the model with the smallest prediction
  error.  The latter is useful for nested models or for
  models with a tuning parameter controlling the complexity
  of the model (e.g., penalized regression).  It selects
  the most parsimonious model whose prediction error is no
  larger than \code{.seFactor} standard errors above the
  prediction error of the best overall model.  Note that
  the models are thereby assumed to be ordered from the
  most parsimonious one to the most complex one.  In
  particular a one-standard-error rule is frequently
  applied.}

  \item{.seFactor}{a numeric value giving a multiplication
  factor of the standard error for the selection of the
  best model.  This is ignored if \code{.selectBest} is
  \code{"min"}.}
}
\value{
  An object of class \code{"perrySelect"} with the
  following components:

  \item{pe}{a data frame containing the estimated
  prediction errors for the models.  In case of more than
  one resampling replication, those are average values over
  all replications.}

  \item{se}{a data frame containing the estimated standard
  errors of the prediction loss for the models.}

  \item{reps}{a data frame containing the estimated
  prediction errors for the models from all replications.
  This is only returned in case of more than one resampling
  replication.}

  \item{splits}{an object giving the data splits used to
  estimate the prediction error of the models.}

  \item{y}{the response.}

  \item{yHat}{a list containing the predicted values for
  the models.  Each list component is again a list
  containing the corresponding predicted values from all
  replications.}

  \item{best}{an integer vector giving the indices of the
  models with the best prediction performance.}

  \item{selectBest}{a character string specifying the
  criterion used for selecting the best model.}

  \item{seFactor}{a numeric value giving the multiplication
  factor of the standard error used for the selection of
  the best model.}
}
\description{
  Combine resampling-based prediction error results for
  various models into one object and select the model with
  the best prediction performance.
}
\details{
  Keep in mind that objects inheriting from class
  \code{"perry"} or \code{"perrySelect"} may contain
  multiple columns of prediction error results.  This is
  the case if the response is univariate but the function
  to compute predictions (usually the
  \code{\link[stats]{predict}} method of the fitted model)
  returns a matrix.

  The \code{.reshape} argument determines how to handle
  such objects.  If \code{.reshape} is \code{FALSE}, all
  objects are required to have the same number of columns
  and the best model for each column is selected.  A
  typical use case for this behavior would be if the
  investigated models contain prediction error results for
  a raw and a reweighted fit.  It might then be of interest
  to researchers to compare the best model for the raw
  estimators with the best model for the reweighted
  estimators.

  If \code{.reshape} is \code{TRUE}, objects with more than
  one column of results are first transformed with
  \code{\link{perryReshape}} to have only one column.  Then
  the best overall model is selected.

  It should also be noted that the argument names of
  \code{.list}, \code{.reshape}, \code{.selectBest} and
  \code{.seFacor} start with a dot to avoid conflicts with
  the argument names used for the objects containing
  prediction error results.
}
\note{
  To ensure comparability, the prediction errors for all
  models are required to be computed from the same data
  splits.
}
\examples{
library("perryExamples")
data("coleman")
set.seed(1234)  # set seed for reproducibility

## set up folds for cross-validation
folds <- cvFolds(nrow(coleman), K = 5, R = 10)

## compare LS, MM and LTS regression

# perform cross-validation for an LS regression model
fitLm <- lm(Y ~ ., data = coleman)
cvLm <- perry(fitLm, splits = folds, 
    cost = rtmspe, trim = 0.1)

# perform cross-validation for an MM regression model
fitLmrob <- lmrob(Y ~ ., data = coleman)
cvLmrob <- perry(fitLmrob, splits = folds, 
    cost = rtmspe, trim = 0.1)

# perform cross-validation for an LTS regression model
fitLts <- ltsReg(Y ~ ., data = coleman)
cvLts <- perry(fitLts, splits = folds, 
    cost = rtmspe, trim = 0.1)

# compare cross-validation results
perrySelect(LS = cvLm, MM = cvLmrob, LTS = cvLts)
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
  \code{\link{perryFit}}, \code{\link{perryTuning}}
}
\keyword{utilities}

