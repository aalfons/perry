# --------------------------------------
# Author: Andreas Alfons
#         Erasmus Universiteit Rotterdam
# --------------------------------------

#' Set up a plot of resampling-based prediction error results
#'
#' Extract and prepare the relevant information for a plot of results of
#' resampling-based prediction error measures.
#'
#' @param object  an object inheriting from class \code{"perry"} or
#' \code{"perrySelect"} that contains prediction error results.
#' @param subset  a character, integer or logical vector indicating the subset
#' of models to be prepared for plotting.
#' @param select  a character, integer or logical vector indicating the columns
#' of prediction error results to be prepared for plotting.
#' @param reps  a logical indicating whether to prepare the results from all
#' replications (\code{TRUE}) or the aggregated results (\code{FALSE}).  The
#' former is suitable for box plots or smooth density plots, while the latter
#' is suitable for dot plots or line plots (see \code{\link{perryPlot}}).
#' @param seFactor  a numeric value giving the multiplication factor of the
#' standard error for displaying error bars in dot plots or line plots.  Error
#' bars in those plots can be suppressed by setting this to \code{NA}.
#' @param \dots  for the \code{"perryTuning"} method, additional arguments to
#' be passed down to the \code{"perrySelect"} method.  For the other methods,
#' additional arguments are currently ignored.
#'
#' @return  An object of class \code{"setupPerryPlot"} with the following
#' components:
#' \describe{
#'   \item{\code{data}}{a data frame containing the following columns:
#'   \describe{
#'     \item{\code{Fit}}{a vector or factor containing the identifiers of the
#'     models.}
#'     \item{\code{Name}}{a factor containing the names of the predictor error
#'     results (not returned in case of only one column of prediction error
#'     results with the default name).}
#'     \item{\code{PE}}{the estimated prediction errors.}
#'     \item{\code{Lower}}{the lower end points of the error bars (only
#'     returned if possible to compute).}
#'     \item{\code{Upper}}{the upper end points of the error bars (only
#'     returned if possible to compute).}
#'   }
#'   }
#'   \item{\code{reps}}{a logical indicating whether the results from all
#'   replications or the aggregated results have been prepared.}
#'   \item{\code{includeSE}}{a logical indicating whether error bars based on
#'   standard errors are available (only returned in case of aggregated
#'   results).}
#'   \item{\code{facets}}{default faceting formula for the plots (not returned
#'   in case of only one column of prediction error results with the default
#'   name).}
#'   \item{\code{tuning}}{a data frame containing the grid of tuning parameter
#'   values for which the prediction error was estimated (only returned for the
#'   \code{"perryTuning"} method).}
#' }
#'
#' @note Duplicate indices in \code{subset} or \code{select} are removed such
#' that all models and prediction error results are unique.
#'
#' @author Andreas Alfons
#'
#' @seealso
#' \code{\link{perryPlot}},
#'
#' \code{\link{perryFit}}, \code{\link{perrySelect}},
#' \code{\link{perryTuning}},
#'
#' \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{autoplot}},
#' \code{\link[graphics]{plot}}
#'
#' @example inst/doc/examples/example-setupPerryPlot.R
#'
#' @keywords utilities
#'
#' @export

setupPerryPlot <- function(object, ...) UseMethod("setupPerryPlot")


#' @rdname setupPerryPlot
#' @method setupPerryPlot perry
#' @export

setupPerryPlot.perry <- function(object, select = NULL, reps = NULL,
                                 seFactor = 1, ...) {
  # initializations
  if (is.null(reps)) reps <- object$splits$R > 1
  else reps <- isTRUE(reps)
  # extract subset of models
  object <- subset(object, select = select)
  if (reps) {
    PE <- object$reps
    if (is.null(PE)) stop("replications not available")
    else PE <- as.data.frame(PE)
  } else PE <- as.data.frame(t(object$pe))
  if (npe(object) == 0) stop("empty prediction error object")
  # stack selected results on top of each other
  fitName <- defaultFitNames(1)
  peName <- defaultNames(1)
  peNames <- peNames(object)
  n <- nrow(PE)
  Fit <- data.frame(Fit = rep.int(fitName, n))
  # no column for conditional plots if there is only one method with default
  # name
  if (isTRUE(peNames == peName)) {
    PE <- cbind(Fit, PE)
    facets <- NULL
  } else {
    PE <- lapply(peNames,
                 function(j) cbind(Fit, Name = rep.int(j, n), PE = PE[, j]))
    PE <- do.call(rbind, PE)
    names(PE) <- c("Fit", "Name", peName)
    facets <- ~ Name
  }
  # add data for error bars unless all replications are requested
  if (reps) includeSE <- NULL
  else {
    includeSE <- !is.null(seFactor) && !is.na(seFactor) &&
      !all(is.na(object$se))
    if (includeSE) {
      halflength <- seFactor * object$se
      PE$Lower <- PE[, peName] - halflength
      PE$Upper <- PE[, peName] + halflength
    }
  }
  # construct object to return
  out <- list(data = PE, reps = reps)
  if (!is.null(includeSE)) out$includeSE <- includeSE
  if (!is.null(facets)) out$facets <- facets
  class(out) <- "setupPerryPlot"
  out
}


#' @rdname setupPerryPlot
#' @method setupPerryPlot perrySelect
#' @export

setupPerryPlot.perrySelect <- function(object, subset = NULL, select = NULL,
                                       reps = NULL, seFactor = NULL, ...) {
  # initializations
  if (is.null(reps)) reps <- object$splits$R > 1
  else reps <- isTRUE(reps)
  if (is.null(seFactor)) seFactor <- object$seFactor
  # extract subset of models
  object <- subset(object, subset = subset, select = select)
  fits <- fits(object)
  if (reps) {
    PE <- object$reps
    if (is.null(PE)) stop("replications not available")
  } else PE <- object$pe
  if (nfits(object) == 0 || npe(object) == 0) {
    stop("empty prediction error object")
  }
  # ensure that models are shown in the correct order and drop unused levels
  # ensure that correct values are shown for a numeric tuning parameter
  if (!is.numeric(PE[, "Fit"])) {
    fits <- fits(object)
    PE$Fit <- factor(PE[, "Fit"], levels = fits)
  }
  # stack selected results on top of each other
  peName <- defaultNames(1)
  peNames <- peNames(object)
  n <- nrow(PE)
  # no column for conditional plots if there is only one column of results
  # with default name
  if (isTRUE(peNames == peName)) facets <- NULL
  else {
    Fit <- PE[, "Fit", drop = FALSE]
    PE <- lapply(peNames,
                 function(j) cbind(Fit, Name = rep.int(j, n), PE = PE[, j]))
    PE <- do.call(rbind, PE)
    names(PE) <- c("Fit", "Name", peName)
    facets <- ~ Name
  }
  # add data for error bars unless all replications are requested
  if (reps) includeSE <- NULL
  else {
    includeSE <- !is.null(seFactor) && !is.na(seFactor) &&
      !all(is.na(object$se[, peNames]))
    if (includeSE) {
      halflength <- seFactor * unlist(object$se[, peNames], use.names = FALSE)
      PE$Lower <- PE[, peName] - halflength
      PE$Upper <- PE[, peName] + halflength
    }
  }
  # construct object to return
  out <- list(data = PE, reps = reps)
  if (!is.null(includeSE)) out$includeSE <- includeSE
  if (!is.null(facets)) out$facets <- facets
  class(out) <- "setupPerryPlot"
  out
}


#' @rdname setupPerryPlot
#' @method setupPerryPlot perryTuning
#' @export

setupPerryPlot.perryTuning <- function(object, ...) {
  #initializations
  tuning <- object$tuning
  # adjust column specifying the model in case of only one tuning parameter
  if (ncol(tuning) == 1) fits(object) <- tuning[, 1]
  # call method for class "perrySelect"
  out <- setupPerryPlot.perrySelect(object, ...)
  # add information on tuning parameters
  out$tuning <- tuning
  out
}
