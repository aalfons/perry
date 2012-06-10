# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' Kernel density plots of resampling-based prediction error results
#' 
#' Produce kernel density plots of results for resampling-based prediction 
#' error measures.
#' 
#' For objects with multiple columns of prediction error results, conditional 
#' kernel density plots are produced.
#' 
#' @method densityplot perry
#' 
#' @param x  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results.
#' @param data  currently ignored.
#' @param subset  a character, integer or logical vector indicating the subset 
#' of models for which to plot the prediction error results.
#' @param select  a character, integer or logical vector indicating the columns 
#' of prediction error results to be plotted.
#' @param \dots  additional arguments to be passed to the \code{"formula"} 
#' method of \code{\link[lattice:histogram]{densityplot}}.
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
#' \code{\link[=bwplot.perry]{bwplot}}, \code{\link[=xyplot.perry]{xyplot}}, 
#' \code{\link[=dotplot.perry]{dotplot}}
#' 
## @example inst/doc/examples/example-densityplot.R
#' 
#' @keywords hplot
#' 
#' @import lattice
#' @export

densityplot.perry <- function(x, data, select = NULL, ...) {
    # initializations
    if(x$splits$R == 1) stop("density plot is not meaningful for one replication")
    # construct data frame in lattice format and call internal function
    PE <- getLatticeData(x, select)
    localDensityplot(PE, ...)
}


#' @rdname densityplot.perry
#' @method densityplot perrySelect
#' @export

densityplot.perrySelect <- function(x, data, subset = NULL, select = NULL, ...) {
    # initializations
    if(x$splits$R == 1) stop("density plot is is not meaningful for one replication")
    # construct data frame in lattice format and call internal function
    PE <- getLatticeData(x, subset, select, numericAsFactor=TRUE)
    localDensityplot(PE, ...)
}


# internal function for density plots
localDensityplot <- function(PE, auto.key = TRUE, xlab = "Prediction error", ...,
        # the following arguments are defined so that they aren't supplied twice
        x, formula, data, groups) {
    # prepare legend
    if(isTRUE(auto.key)) {
        auto.key <- list(points=TRUE, lines=TRUE)
    } else if(is.list(auto.key)) {
        if(is.null(auto.key$points)) auto.key$points <- TRUE
        if(is.null(auto.key$lines)) auto.key$lines <- TRUE
    }
    # construct formula for call to densityplot()
    peNames <- names(PE)
    conditional <- if("Name" %in% peNames) "Name" else NULL
    f <- getFormula("", "PE", conditional)
    # call densityplot()
    if("Fit" %in% peNames) {
        # this is ugly, but avoids NOTE in R CMD CHECK
        command <- paste("densityplot(f, data=PE, groups=Fit,", 
            "auto.key=auto.key, xlab=xlab, ...)")
        eval(parse(text=command))
    } else {
        # no NOTE in this case
        densityplot(f, data=PE, auto.key=auto.key, xlab=xlab, ...)
    }
}
