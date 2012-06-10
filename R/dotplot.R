# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @S3method dotplot perry
dotplot.perry <- function(x, data, select = NULL, seFactor = NA, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, select, reps=FALSE, seFactor=seFactor)
    localDotplot(tmp$PE, tmp$lower, tmp$upper, ...)
}

#' @S3method dotplot perrySelect
dotplot.perrySelect <- function(x, data, subset = NULL, select = NULL, 
        seFactor = x$seFactor, ...) {
    # construct data frame in lattice format and call internal function
    tmp <- getLatticeData(x, subset, select, reps=FALSE, 
        seFactor=seFactor, numericAsFactor=TRUE)
    localDotplot(tmp$PE, tmp$lower, tmp$upper, ...)
}


# internal function for dot plots
localDotplot <- function(PE, lower, upper, horizontal = TRUE, 
        xlab = NULL, ylab = NULL, ...,
        # the following arguments are defined so that they aren't supplied twice
        x, formula, data, groups) {
    # construct formula for call to xyplot()
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
    # call stripplot() since dotplot() does not pass the 'subscripts' to the 
    # prepanel function
    stripplot(f, data=PE, lower=lower, upper=upper, horizontal=horizontal, 
        prepanel=prepanelDotplot, panel=panelDotplot, xlab=xlab, ylab=ylab, 
        ...)
}

# prepanel function
prepanelDotplot <- function(x, y, lower, upper, horizontal, subscripts, ...) {
    tmp <- c(lower[subscripts], if(horizontal) x else y, upper[subscripts])
    tmp <- tmp[is.finite(tmp)]
    if(length(tmp) > 0) {
        lim <- range(tmp, finite=TRUE)
        if(horizontal) list(xlim=lim) else list(ylim=lim)
    } else list()
}

# panel function
panelDotplot <- function(x, y, lower, upper, horizontal, subscripts, 
        col = trellis.par.get("dot.symbol")$col, angle=90, 
        length=0.5, unit="lines", ends, type, lty, lwd, ...) {
    # initializations
    dot.line <- trellis.par.get("dot.line")
    box.umbrella <- trellis.par.get("box.umbrella")
    if(missing(lty) || length(lty) == 0) {
        lty <- c(dot.line$lty, box.umbrella$lty)
    } else if(length(lty == 1)) lty = c(lty, box.umbrella$lty)
    if(missing(lwd) || length(lwd) == 0) {
        lwd <- c(dot.line$lwd, box.umbrella$lwd)
    } else if(length(lwd == 1)) lwd = c(lwd, box.umbrella$lwd)
    # create plot
    panel.dotplot(x, y, horizontal=horizontal, subscripts=subscripts, 
        col=col, lty=lty[1], lwd=lwd[1], ...)
    if(horizontal) {
        panel.arrows(lower[subscripts], y, upper[subscripts], y, 
            angle=angle, length=length, unit=unit, ends="both", 
            col=col, lty=lty[2], lwd=lwd[2], ...)
    } else {
        panel.arrows(x, lower[subscripts], x, upper[subscripts], 
            angle=angle, length=length, unit=unit, ends="both", 
            col=col, lty=lty[2], lwd=lwd[2], ...)
    }
}
