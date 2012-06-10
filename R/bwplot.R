# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @S3method bwplot perry
bwplot.perry <- function(x, data, select = NULL, ...) {
    # initializations
    if(x$splits$R == 1) stop("box plot is not meaningful for one replication")
    # construct data frame in lattice format and call internal function
    PE <- getLatticeData(x, select)
    localBwplot(PE, ...)
}

#' @S3method bwplot perrySelect
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
