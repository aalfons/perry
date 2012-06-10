\name{plot.perry}
\alias{plot.perry}
\alias{plot.perrySelect}
\title{Plot resampling-based prediction error results}
\usage{
  \method{plot}{perry} (x,
    method = c("bwplot", "densityplot", "xyplot", "dotplot"),
    ...)

  \method{plot}{perrySelect} (x,
    method = c("bwplot", "densityplot", "xyplot", "dotplot"),
    ...)
}
\arguments{
  \item{x}{an object inheriting from class \code{"perry"}
  or \code{"perrySelect"} that contains prediction error
  results.}

  \item{method}{a character string specifying the type of
  plot.  Possible values are \code{"bwplot"} to create a
  box-and-whisker plot via
  \code{\link[=bwplot.perry]{bwplot}} (the default),
  \code{"densityplot"} to create a kernel density plot via
  \code{\link[=densityplot.perry]{densityplot}},
  \code{"xyplot"} to plot the (average) results for each
  model via \code{\link[=xyplot.perry]{xyplot}}, or
  \code{"dotplot"} to create a similar dot plot via
  \code{\link[=dotplot.perry]{dotplot}}.  Note that the
  first two plots are only meaningful in case of repeated
  resampling.  The default is to use \code{"bwplot"} in
  case of repeated resampling and \code{"xyplot"}
  otherwise.  In any case, partial string matching allows
  supply abbreviations of the accepted values.}

  \item{\dots}{additional arguments to be passed down.}
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
  Plot results of resampling-based prediction error
  measures.
}
\details{
  For objects with multiple columns of prediction error
  results, conditional plots are produced.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perrySelect}},
  \code{\link{perryTuning}},
  \code{\link[=bwplot.perry]{bwplot}},
  \code{\link[=densityplot.perry]{densityplot}},
  \code{\link[=xyplot.perry]{xyplot}},
  \code{\link[=dotplot.perry]{dotplot}}
}
\keyword{hplot}

