# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' X-Y plots of resampling-based prediction error results
#' 
#' Plot the (average) results for resampling-based prediction error measures on 
#' the \eqn{y}-axis against the respective models on the \eqn{x}-axis.
#' 
#' For objects with multiple columns of prediction error results, conditional 
#' plots are produced.
#' 
#' In most situations, the default behavior is to represent the 
#' prediction error results for each model by a vertical line segment (i.e., to 
#' call the default method of \code{\link[lattice:xyplot]{xyplot}} with 
#' \code{type = "h"}).  However, the behavior is different for objects of class 
#' \code{"perryTuning"} with only one numeric tuning parameter.  In that 
#' situation, the prediction error results are plotted against the values of 
#' the tuning parameter as a connected line (i.e., by using \code{type = "b"}).
#' 
#' The default behavior can of course be overridden by supplying the 
#' \code{type} argument (a full list of accepted values can be found in the 
#' help file of \code{\link[lattice:panel.xyplot]{panel.xyplot}}).
#' 
#' @method xyplot perry
#' 
#' @param x  an object inheriting from class \code{"perry"} or 
#' \code{"perrySelect"} that contains prediction error results.
#' @param data  currently ignored.
#' @param subset  a character, integer or logical vector indicating the subset 
#' of models for which to plot the prediction error results.
#' @param select  a character, integer or logical vector indicating the columns 
#' of prediction error results to be plotted.
#' @param seFactor  a numeric value giving the multiplication factor of the 
#' standard error for displaying error bars.  Error bars can be suppressed by 
#' setting this to \code{NA}.
#' @param \dots  additional arguments to be passed to the \code{"formula"} 
#' method of \code{\link[lattice:xyplot]{xyplot}}.
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
#' \code{\link[=dotplot.perry]{dotplot}}, \code{\link[=bwplot.perry]{bwplot}}, 
#' \code{\link[=densityplot.perry]{densityplot}}
#' 
## @example inst/doc/examples/example-xyplot.R
#' 
#' @keywords hplot
#' 
#' @import lattice
#' @export

xyplot.perry <- function(x, data, select = NULL, seFactor = NA, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, select, reps=FALSE, seFactor=seFactor)
    localXyplot(tmp$PE, tmp$lower, tmp$upper, ...)
}


#' @rdname xyplot.perry
#' @method xyplot perrySelect
#' @export

xyplot.perrySelect <- function(x, data, subset = NULL, select = NULL, 
        seFactor = x$seFactor, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, subset, select, reps=FALSE, 
        seFactor=seFactor, numericAsFactor=TRUE)
    localXyplot(tmp$PE, tmp$lower, tmp$upper, ...)
}


#' @rdname xyplot.perry
#' @method xyplot perryTuning
#' @export

xyplot.perryTuning <- function(x, data, subset = NULL, select = NULL, 
        seFactor = x$seFactor, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, subset, select, reps=FALSE, seFactor=seFactor)
    localXyplot(tmp$PE, tmp$lower, tmp$upper, x$tuning, ...)
}


# internal function for x-y plots
localXyplot <- function(PE, lower, upper, tuning = NULL, type, 
        xlab, ylab = "Prediction error", ...,
        # the following arguments are defined so that they aren't supplied twice
        x, formula, data, groups) {
    # construct formula for call to xyplot()
    peNames <- names(PE)
    if(!("Fit" %in% peNames)) {
        PE$Fit <- factor(rep.int(defaultFitNames(1), nrow(PE)))
    }
    conditional <- if("Name" %in% peNames) "Name" else NULL
    f <- getFormula("PE", "Fit", conditional)
    # default plot type x-axis label
    if(!is.null(tuning) && length(names(tuning)) == 1 && is.numeric(PE$Fit)) {
        if(missing(type)) type <- "b"
        if(missing(xlab)) xlab <- names(tuning)
    } else {
        if(missing(type)) type <- c("h", "p")
        if(missing(xlab)) xlab <- NULL
    }
    # call xyplot()
    xyplot(f, data=PE, lower=lower, upper=upper, prepanel=prepanelXyplot, 
        panel=panelXyplot, type=type, xlab=xlab, ylab=ylab, ...)
}

# prepanel function
prepanelXyplot <- function(x, y, lower, upper, subscripts, ...) {
    tmp <- c(lower[subscripts], y, upper[subscripts])
    tmp <- tmp[is.finite(tmp)]
    if(length(tmp) > 0) {
        lim <- range(tmp, finite=TRUE)
        list(ylim=lim)
    } else list()
}

# panel function
panelXyplot <- function(x, y, lower, upper, subscripts, 
        col = plot.line$col, angle=90, length=0.5, 
        unit="lines", ends, type, lty, lwd, ...) {
    # initializations
    plot.line <- trellis.par.get("plot.line")
    box.umbrella <- trellis.par.get("box.umbrella")
    if(missing(lty) || length(lty) == 0) {
        lty <- c(plot.line$lty, box.umbrella$lty)
    } else if(length(lty == 1)) lty = c(lty, box.umbrella$lty)
    if(missing(lwd) || length(lwd) == 0) {
        lwd <- c(plot.line$lwd, box.umbrella$lwd)
    } else if(length(lwd == 1)) lwd = c(lwd, box.umbrella$lwd)
    # create plot
    panel.xyplot(x, y, subscripts=subscripts, type=type, col=col, 
        lty=lty[1], lwd=lwd[1], ...)
    panel.arrows(x, lower[subscripts], x, upper[subscripts], 
        angle=angle, length=length, unit=unit, ends="both", 
        col=col, lty=lty[2], lwd=lwd[2], ...)
}
