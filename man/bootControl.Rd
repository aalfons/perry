\name{bootControl}
\alias{bootControl}
\title{Control object for bootstrap samples}
\usage{
  bootControl(R = 1, type = c("0.632", "out-of-bag"),
    grouping = NULL)
}
\arguments{
  \item{R}{an integer giving the number of bootstrap
  samples.}

  \item{type}{a character string specifying a bootstrap
  estimator.  Possible values are \code{"0.632"} (the
  default), or \code{"out-of-bag"}.}

  \item{grouping}{a factor specifying groups of
  observations.}
}
\value{
  An object of class \code{"bootSamples"} with the
  following components:

  \item{R}{an integer giving the number of bootstrap
  samples.}

  \item{type}{a character string specifying the type of
  bootstrap estimator.}

  \item{grouping}{if supplied, a factor specifying groups
  of observations.  The groups will then be resampled
  rather than individual observations such that all
  observations within a group belong either to the
  bootstrap sample or the test data.}
}
\description{
  Generate an object that controls how to draw bootstrap
  samples and which bootstrap estimator of prediction error
  to compute.
}
\examples{
set.seed(1234)  # set seed for reproducibility
perrySplits(20, bootControl())
perrySplits(20, bootControl(R = 10))
}
\author{
  Andreas Alfons
}
\references{
  Efron, B. (1983) Estimating the error rate of a
  prediction rule: improvement on cross-validation.
  \emph{Journal of the American Statistical Association},
  \bold{78}(382), 316--331.
}
\seealso{
  \code{\link{perrySplits}}, \code{\link{bootSamples}},
  \code{\link{foldControl}}, \code{\link{splitControl}}
}
\keyword{utilities}

