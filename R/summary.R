# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Summarize prediction error results
#' 
#' Produce a summary of results from (repeated) \eqn{K}-fold cross-validation, 
#' (repeated) random splitting (also known as random subsampling or Monte Carlo 
#' cross-validation), or the bootstrap.  
#' 
#' @method summary perry
#' 
#' @param object  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results (note that the 
#' latter includes objects of class \code{"perryTuning"}).
#' @param \dots  currently ignored.
#' 
#' @return 
#' An object of class \code{"summary.perry"}, \code{"summary.perrySelect"} or 
#' \code{"summary.perryTuning"}, depending on the class of \code{object}.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perryFit}}, \code{\link{perrySelect}}, 
#' \code{\link{perryTuning}}, \code{\link{summary}}
#' 
## @example inst/doc/examples/example-summary.R
#' 
#' @keywords utilities
#' 
#' @export

summary.perry <- function(object, ...) {
    pe <- aggregate(object, summary)
    out <- list(splits=object$splits, pe=pe)
    class(out) <- "summary.perry"
    out
}


#' @rdname summary.perry
#' @method summary perrySelect
#' @export

summary.perrySelect <- function(object, ...) {
    pe <- aggregate(object, summary)
    out <- list(splits=object$splits, best=object$best, pe=pe)
    class(out) <- "summary.perrySelect"
    out
}


#' @rdname summary.perry
#' @method summary perryTuning
#' @export

summary.perryTuning <- function(object, ...) {
    out <- summary.perrySelect(object, ...)
    out <- list(splits=object$splits, tuning=object$tuning, 
        best=out$best, pe=out$pe)
    class(out) <- c("summary.perryTuning", class(out))
    out
}
