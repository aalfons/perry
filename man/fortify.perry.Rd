\name{fortify.perry}
\alias{fortify.perry}
\alias{fortify.perrySelect}
\alias{fortify.perryTuning}
\title{Convert resampling-based prediction error results into a data frame for
plotting}
\usage{
  \method{fortify}{perry} (model, data, select = NULL,
    reps = model$splits$R > 1, seFactor = NA, ...)

  \method{fortify}{perrySelect} (model, data,
    subset = NULL, select = NULL,
    reps = model$splits$R > 1, seFactor = model$seFactor,
    ...)

  \method{fortify}{perryTuning} (model, data, ...)
}
\arguments{
  \item{model}{an object inheriting from class
  \code{"perry"} or \code{"perrySelect"} that contains
  prediction error results.}

  \item{data}{currently ignored.}

  \item{subset}{a character, integer or logical vector
  indicating the subset of models to be converted.}

  \item{select}{a character, integer or logical vector
  indicating the columns of prediction error results to be
  converted.}

  \item{reps}{a logical indicating whether to convert the
  results from all replications (\code{TRUE}) or the
  aggregated results (\code{FALSE}).  The former is
  suitable for box plots or smooth density plots, while the
  latter is suitable for dot plots or line plots (see
  \code{\link{perryPlot}}).}

  \item{seFactor}{a numeric value giving the multiplication
  factor of the standard error for displaying error bars in
  dot plots or line plots.  Error bars in those plots can
  be suppressed by setting this to \code{NA}.}

  \item{\dots}{for the \code{"perryTuning"} method,
  additional arguments to be passed down to the
  \code{"perrySelect"} method.  For the other methods,
  additional arguments are currently ignored.}
}
\value{
  A data frame containing the columns listed below, as well
  as additional information stored in the attribute
  \code{"facets"} (default faceting formula for the plots).

  \item{Fit}{a vector or factor containing the identifiers
  of the models.}

  \item{Name}{a factor containing the names of the
  predictor error results (not returned in case of only one
  column of prediction error results with the default
  name).}

  \item{PE}{the estimated prediction errors.}

  \item{Lower}{the lower end points of the error bars (only
  returned if \code{reps} is \code{FALSE}).}

  \item{Upper}{the upper end points of the error bars (only
  returned if \code{reps} is \code{FALSE}).}
}
\description{
  Extract all necessary information for plotting from
  resampling-based prediction error results and store it in
  a data frame.
}
\note{
  Duplicate indices in \code{subset} or \code{select} are
  removed such that all models and prediction error results
  are unique.
}
\author{
  Andreas Alfons
}
\seealso{
  \code{\link[ggplot2]{fortify}}, \code{\link{perryPlot}},
  \code{\link{perryFit}}, \code{\link{perrySelect}},
  \code{\link{perryTuning}}
}
\keyword{utilities}

