\name{summary.perry}
\alias{summary.perry}
\alias{summary.perrySelect}
\alias{summary.perryTuning}
\title{Summarize prediction error results}
\usage{
  \method{summary}{perry} (object, ...)

  \method{summary}{perrySelect} (object, ...)

  \method{summary}{perryTuning} (object, ...)
}
\arguments{
  \item{object}{an object inheriting from class
  \code{"perry"} or \code{"perrySelect"} that contains
  prediction error results (note that the latter includes
  objects of class \code{"perryTuning"}).}

  \item{\dots}{currently ignored.}
}
\value{
  An object of class \code{"summary.perry"},
  \code{"summary.perrySelect"} or
  \code{"summary.perryTuning"}, depending on the class of
  \code{object}.
}
\description{
  Produce a summary of results from (repeated) \eqn{K}-fold
  cross-validation, Monte Carlo cross-validation(also known
  as random splitting or random subsampling), or the
  bootstrap.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perrySelect}},
  \code{\link{perryTuning}}, \code{\link{summary}}
}
\keyword{utilities}

