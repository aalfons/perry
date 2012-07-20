\name{bootPE}
\alias{bootPE}
\title{Bootstrap prediction error estimation for fitted models}
\usage{
  bootPE(object, R = 1,
    bootType = c("0.632", "out-of-bag"), grouping = NULL,
    samples = NULL, ...)
}
\arguments{
  \item{object}{the fitted model for which to estimate the
  prediction error.}

  \item{R}{an integer giving the number of bootstrap
  samples.}

  \item{bootType}{a character string specifying a bootstrap
  estimator.  Possible values are \code{"0.632"} (the
  default), or \code{"out-of-bag"}.}

  \item{grouping}{a factor specifying groups of
  observations.  If supplied, the groups are resampled
  rather than individual observations such that all
  observations within a group belong either to the
  bootstrap sample or the test data.}

  \item{samples}{an object of class \code{"bootSamples"}
  (as returned by \code{\link{bootSamples}}) or a control
  object of class \code{"bootControl"} (see
  \code{\link{bootControl}}) defining bootstrap samples.
  If supplied, this is preferred over the arguments for
  generating bootstrap samples.}

  \item{\dots}{additional arguments to be passed down to
  \code{\link{perry}}.}
}
\value{
  An object of class \code{"perry"} as returned by
  \code{\link{perry}}.
}
\description{
  Estimate the prediction error of a fitted model via the
  bootstrap.  This works for any model for which a
  \code{\link{perry}} method is available.
}
\examples{
# load data and fit an LS regression model
data("mtcars")
fit <- lm(mpg ~ wt + cyl, data=mtcars)

# perform bootstrap prediction error estimation
bootPE(fit, R = 10, bootType = "0.632", seed = 1234)
bootPE(fit, R = 10, bootType = "out-of-bag", seed = 1234)
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perry}}, \code{\link{repCV}},
  \code{\link{repRS}}
}
\keyword{utilities}

