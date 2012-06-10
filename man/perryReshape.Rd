\name{perryReshape}
\alias{perryReshape}
\alias{perryReshape.perry}
\alias{perryReshape.perrySelect}
\title{Reshape resampling-based prediction error results}
\usage{
  perryReshape(x, ...)

  \method{perryReshape}{perry} (x,
    selectBest = c("min", "hastie"), seFactor = 1, ...)

  \method{perryReshape}{perrySelect} (x,
    selectBest = c("min", "hastie"), seFactor = 1, ...)
}
\arguments{
  \item{x}{an object inheriting from class \code{"perry"}
  or \code{"perrySelect"} that contains prediction error
  results.}

  \item{selectBest}{a character string specifying a
  criterion for selecting the best model.  Possible values
  are \code{"min"} (the default) or \code{"hastie"}.  The
  former selects the model with the smallest prediction
  error.  The latter is useful for nested models or for
  models with a tuning parameter controlling the complexity
  of the model (e.g., penalized regression).  It selects
  the most parsimonious model whose prediction error is no
  larger than \code{seFactor} standard errors above the
  prediction error of the best overall model.  Note that
  the models are thereby assumed to be ordered from the
  most parsimonious one to the most complex one.  In
  particular a one-standard-error rule is frequently
  applied.}

  \item{seFactor}{a numeric value giving a multiplication
  factor of the standard error for the selection of the
  best model.  This is ignored if \code{selectBest} is
  \code{"min"}.}

  \item{\dots}{additional arguments to be passed down.}
}
\value{
  An object of class \code{"perrySelect"} with the
  following components:

  \item{splits}{an object giving the data splits used to
  estimate the prediction error.}

  \item{best}{an integer giving the index of the model with
  the best prediction performance.}

  \item{pe}{a data frame containing the estimated
  prediction errors for the models.  In case of more than
  one resampling replication, those are average values over
  all replications.}

  \item{se}{a data frame containing the estimated standard
  errors of the prediction loss for the models.}

  \item{selectBest}{a character string specifying the
  criterion used for selecting the best model.}

  \item{seFactor}{a numeric value giving the multiplication
  factor of the standard error used for the selection of
  the best model.}

  \item{reps}{a data frame containing the estimated
  prediction errors for the models from all replications.
  This is only returned in case of more than one resampling
  replication.}
}
\description{
  Reshape resampling-based prediction error results into an
  object of class \code{"perrySelect"} with only one column
  of results.
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
  \code{\link{perryTuning}}
}
\keyword{utilities}

