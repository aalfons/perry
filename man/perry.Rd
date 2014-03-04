\name{perry}
\alias{perry}
\title{Resampling-based prediction error for fitted models}
\usage{
  perry(object, ...)
}
\arguments{
  \item{object}{the fitted model for which to estimate the
  prediction error.}

  \item{\dots}{additional arguments to be passed down to
  methods.}
}
\description{
  Generic function to estimate the prediction error of a
  fitted model via (repeated) \eqn{K}-fold
  cross-validation, (repeated) random splitting (also known
  as random subsampling or Monte Carlo cross-validation),
  or the bootstrap.
}
\details{
  The idea is that developers write easy-to-use methods for
  end users to leverage the prediction error estimation
  framework for their models.  A typical \code{perry}
  method consists of the following two parts: first the
  data are extracted from the model, then function
  \code{\link{perryFit}} is called to perform prediction
  error estimation.  The programming effort of implementing
  prediction error estimation for a certain model is thus
  greatly reduced.

  Examples for methods are available in package
  perryExamples (see
  \code{\link[perryExamples]{perry-methods}}).
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}},
  \code{\link[perryExamples]{perry-methods}}
}
\keyword{utilities}

