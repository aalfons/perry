# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Box-and-whisker plots of resampling-based prediction error results
#' 
#' Produce box-and-whisker plots of results for resampling-based prediction 
#' error measures.
#' 
#' For objects with multiple columns of prediction error results, conditional 
#' box-and-whisker plots are produced.
#' 
#' @method bwplot perry
#' 
#' @param x  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results.
#' @param data  currently ignored.
#' @param subset  a character, integer or logical vector indicating the subset 
#' of models for which to plot the prediction error results.
#' @param select  a character, integer or logical vector indicating the columns 
#' of prediction error results to be plotted.
#' @param \dots  additional arguments to be passed to the \code{"formula"} 
#' method of \code{\link[lattice:xyplot]{bwplot}}.
#' 
#' @return An object of class \code{"trellis"} is returned invisibly.  The 
#' \code{\link[lattice:update.trellis]{update}} method can be used to update 
#' components of the object and the \code{\link[lattice:print.trellis]{print}} 
#' method (usually called by default) will plot it on an appropriate plotting 
#' device.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{perryFit}}, \code{\link{perrySelect}}, 
#' \code{\link{perryTuning}}, \code{\link[=plot.perry]{plot}}, 
#' \code{\link[=densityplot.perry]{densityplot}}, 
#' \code{\link[=xyplot.perry]{xyplot}}, \code{\link[=dotplot.perry]{dotplot}}
#' 
## @example inst/doc/examples/example-bwplot.R
#' 
#' @keywords hplot
#' 
#' @import lattice
#' @export

bwplot.perry <- function(x, data, select = NULL, ...) {
    # initializations
    if(x$splits$R == 1) stop("box plot is not meaningful for one replication")
    # construct data frame in lattice format and call internal function
    PE <- getLatticeData(x, select)
    localBwplot(PE, ...)
}


#' @rdname bwplot.perry
#' @method bwplot perrySelect
#' @export

bwplot.perrySelect <- function(x, data, subset = NULL, select = NULL, ...) {
    # initializations
    if(x$splits$R == 1) stop("box plot is not meaningful for one replication")
    # construct data frame in lattice format and call internal function
    PE <- getLatticeData(x, subset, select, numericAsFactor=TRUE)
    localBwplot(PE, ...)
}


## internal function for box plots
localBwplot <- function(PE, horizontal = TRUE, xlab = NULL, ylab = NULL, ..., 
        # the following arguments are defined so that they aren't supplied twice
        x, formula, data, groups) {
    # construct formula for call to bwplot()
    peNames <- names(PE)
    haveFit <- "Fit" %in% peNames
    if(horizontal) {
        left <- if(haveFit) "Fit" else ""
        right <- "PE"
        if(missing(xlab)) xlab <- "Prediction error"
    } else {
        if(!haveFit) PE$Fit <- rep.int(defaultFitNames(1), nrow(PE))
        left <- "PE"
        right <- "Fit"
        if(missing(ylab)) ylab <- "Prediction error"
    }
    conditional <- if("Name" %in% peNames) "Name" else NULL
    f <- getFormula(left, right, conditional)
    # call bwplot()
    bwplot(f, data=PE, horizontal=horizontal, xlab=xlab, ylab=ylab, ...)
}
