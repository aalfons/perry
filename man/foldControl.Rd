\name{foldControl}
\alias{foldControl}
\title{Control object for cross-validation folds}
\usage{
  foldControl(K = 5, R = 1,
    type = c("random", "consecutive", "interleaved"),
    grouping = NULL)
}
\arguments{
  \item{K}{an integer giving the number of folds into which
  the observations should be split (the default is five).}

  \item{R}{an integer giving the number of replications for
  repeated \eqn{K}-fold cross-validation.}

  \item{type}{a character string specifying the type of
  folds to be generated.  Possible values are
  \code{"random"} (the default), \code{"consecutive"} or
  \code{"interleaved"}.}

  \item{grouping}{a factor specifying groups of
  observations.}
}
\value{
  An object of class \code{"foldControl"} with the
  following components:

  \item{K}{an integer giving the number of folds.  A value
  of \code{K} equal to the number of observations or groups
  yields leave-one-out cross-validation.}

  \item{R}{an integer giving the number of replications.
  This will be ignored for for leave-one-out
  cross-validation and other non-random splits of the
  data.}

  \item{type}{a character string specifying the type of
  folds.}

  \item{grouping}{if supplied, a factor specifying groups
  of observations.  The data will then be split according
  to the groups rather than individual observations such
  that all observations within a group belong to the same
  fold.}
}
\description{
  Generate an object that controls how to split \eqn{n}
  observations or groups of observations into \eqn{K} folds
  to be used for (repeated) \eqn{K}-fold cross-validation.
  \eqn{K} should thereby be chosen such that all folds are
  of approximately equal size.
}
\examples{
set.seed(1234)  # set seed for reproducibility
perrySplits(20, foldControl(K = 5))
perrySplits(20, foldControl(K = 5, R = 10))
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perrySplits}}, \code{\link{splitControl}},
  \code{\link{cvFolds}}
}
\keyword{utilities}

