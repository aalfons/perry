\name{dotplot.perry}
\alias{dotplot.perry}
\alias{dotplot.perrySelect}
\title{Dot plots of resampling-based prediction error results}
\usage{
  \method{dotplot}{perry} (x, data, select = NULL,
    seFactor = NA, ...)

  \method{dotplot}{perrySelect} (x, data, subset = NULL,
    select = NULL, seFactor = x$seFactor, ...)
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

  \item{seFactor}{a numeric value giving the multiplication
  factor of the standard error for displaying error bars.
  Error bars can be suppressed by setting this to
  \code{NA}.}

  \item{\dots}{additional arguments to be passed to the
  \code{"formula"} method of
  \code{\link[lattice:xyplot]{dotplot}}.}
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
  Produce dot plots of (average) results for
  resampling-based prediction error measures.
}
\details{
  For objects with multiple columns of prediction error
  results, conditional dot plots are produced.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perrySelect}},
  \code{\link{perryTuning}},
  \code{\link[=plot.perry]{plot}},
  \code{\link[=xyplot.perry]{xyplot}},
  \code{\link[=bwplot.perry]{bwplot}},
  \code{\link[=densityplot.perry]{densityplot}}
}
\keyword{hplot}

