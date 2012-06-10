# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @S3method xyplot perry
xyplot.perry <- function(x, data, select = NULL, seFactor = NA, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, select, reps=FALSE, seFactor=seFactor)
    localXyplot(tmp$PE, tmp$lower, tmp$upper, ...)
}

#' @S3method xyplot perrySelect
xyplot.perrySelect <- function(x, data, subset = NULL, select = NULL, 
        seFactor = x$seFactor, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, subset, select, reps=FALSE, 
        seFactor=seFactor, numericAsFactor=TRUE)
    localXyplot(tmp$PE, tmp$lower, tmp$upper, ...)
}


#' @S3method xyplot perryTuning
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
