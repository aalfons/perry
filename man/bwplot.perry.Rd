\name{bwplot.perry}
\alias{bwplot.perry}
\alias{bwplot.perrySelect}
\title{Box-and-whisker plots of resampling-based prediction error results}
\usage{
  \method{bwplot}{perry} (x, data, select = NULL, ...)

  \method{bwplot}{perrySelect} (x, data, subset = NULL,
    select = NULL, ...)
}
\arguments{
  \item{x}{an object inheriting from class \code{"perry"}
  or \code{"perrySelect"} that contains prediction error
  results.}

  \item{data}{currently ignored.}

  \item{subset}{a character, integer or logical vector
  indicating the subset of models for which to plot the
  prediction error results.}

  \item{select}{a character, integer or logical vector
  indicating the columns of prediction error results to be
  plotted.}

  \item{\dots}{additional arguments to be passed to the
  \code{"formula"} method of
  \code{\link[lattice:xyplot]{bwplot}}.}
}
\value{
  An object of class \code{"trellis"} is returned
  invisibly.  The
  \code{\link[lattice:update.trellis]{update}} method can
  be used to update components of the object and the
  \code{\link[lattice:print.trellis]{print}} method
  (usually called by default) will plot it on an
  appropriate plotting device.
}
\description{
  Produce box-and-whisker plots of results for
  resampling-based prediction error measures.
}
\details{
  For objects with multiple columns of prediction error
  results, conditional box-and-whisker plots are produced.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perrySelect}},
  \code{\link{perryTuning}},
  \code{\link[=plot.perry]{plot}},
  \code{\link[=densityplot.perry]{densityplot}},
  \code{\link[=xyplot.perry]{xyplot}},
  \code{\link[=dotplot.perry]{dotplot}}
}
\keyword{hplot}

