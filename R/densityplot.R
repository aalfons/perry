# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @S3method densityplot perry
densityplot.perry <- function(x, data, select = NULL, ...) {
    # initializations
    if(x$splits$R == 1) stop("density plot is not meaningful for one replication")
    # construct data frame in lattice format and call internal function
    PE <- getLatticeData(x, select)
    localDensityplot(PE, ...)
}

#' @S3method densityplot perrySelect
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
