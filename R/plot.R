# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## workhorse function
plotPE <- function(x, 
        method = c("bwplot", "densityplot", "xyplot", "dotplot"), ...) {
    ## initializations
    if(x$splits$R == 1) {
        choices <- eval(formals(sys.function())[["method"]])
        if(identical(method, choices)) {
            method <- "xyplot"
        } else method <- match.arg(method, c("xyplot", "dotplot"))
    } else method <- match.arg(method)
    ## call plot function
    if(method == "bwplot") {
        bwplot(x, ...)
    } else if(method == "densityplot") {
        densityplot(x, ...)
    } else if(method == "xyplot") {
        xyplot(x, ...)
    } else dotplot(x, ...)
}


#' Plot resampling-based prediction error results
#' 
#' Plot results of resampling-based prediction error measures.
#' 
#' For objects with multiple columns of prediction error results, conditional 
#' plots are produced.
#' 
#' @method plot perry
#' 
#' @param x  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results.
#' @param method  a character string specifying the type of plot.  Possible 
#' values are \code{"bwplot"} to create a box-and-whisker plot via 
#' \code{\link[=bwplot.perry]{bwplot}} (the default), \code{"densityplot"} to 
#' create a kernel density plot via 
#' \code{\link[=densityplot.perry]{densityplot}}, \code{"xyplot"} to plot the 
#' (average) results for each model via \code{\link[=xyplot.perry]{xyplot}}, or 
#' \code{"dotplot"} to create a similar dot plot via 
#' \code{\link[=dotplot.perry]{dotplot}}.  Note that the first two plots are 
#' only meaningful in case of repeated resampling.  The default is to use 
#' \code{"bwplot"} in case of repeated resampling and \code{"xyplot"} 
#' otherwise.  In any case, partial string matching allows supply abbreviations 
#' of the accepted values.
#' @param \dots  additional arguments to be passed down.
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
#' \code{\link{perryTuning}}, \code{\link[=bwplot.perry]{bwplot}}, 
#' \code{\link[=densityplot.perry]{densityplot}}, 
#' \code{\link[=xyplot.perry]{xyplot}}, \code{\link[=dotplot.perry]{dotplot}}
#' 
## @example inst/doc/examples/example-plot.R
#' 
#' @keywords hplot
#' 
#' @export

plot.perry <- plotPE


#' @rdname plot.perry
#' @method plot perrySelect
#' @export

plot.perrySelect <- plotPE
