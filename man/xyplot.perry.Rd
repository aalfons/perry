\name{xyplot.perry}
\alias{xyplot.perry}
\alias{xyplot.perrySelect}
\alias{xyplot.perryTuning}
\title{X-Y plots of resampling-based prediction error results}
\usage{
  \method{xyplot}{perry} (x, data, select = NULL,
    seFactor = NA, ...)

  \method{xyplot}{perrySelect} (x, data, subset = NULL,
    select = NULL, seFactor = x$seFactor, ...)

  \method{xyplot}{perryTuning} (x, data, subset = NULL,
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
  \code{\link[lattice:xyplot]{xyplot}}.}
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
  Plot the (average) results for resampling-based
  prediction error measures on the \eqn{y}-axis against the
  respective models on the \eqn{x}-axis.
}
\details{
  For objects with multiple columns of prediction error
  results, conditional plots are produced.

  In most situations, the default behavior is to represent
  the prediction error results for each model by a vertical
  line segment (i.e., to call the default method of
  \code{\link[lattice:xyplot]{xyplot}} with \code{type =
  "h"}).  However, the behavior is different for objects of
  class \code{"perryTuning"} with only one numeric tuning
  parameter.  In that situation, the prediction error
  results are plotted against the values of the tuning
  parameter as a connected line (i.e., by using \code{type
  = "b"}).

  The default behavior can of course be overridden by
  supplying the \code{type} argument (a full list of
  accepted values can be found in the help file of
  \code{\link[lattice:panel.xyplot]{panel.xyplot}}).
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link{perryFit}}, \code{\link{perrySelect}},
  \code{\link{perryTuning}},
  \code{\link[=plot.perry]{plot}},
  \code{\link[=dotplot.perry]{dotplot}},
  \code{\link[=bwplot.perry]{bwplot}},
  \code{\link[=densityplot.perry]{densityplot}}
}
\keyword{hplot}

